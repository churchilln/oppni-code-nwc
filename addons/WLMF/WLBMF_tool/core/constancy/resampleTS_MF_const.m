function [coefB]=resampleTS_MF_const(coef, nj, B1, B2, Block, J1flag);
% function [coefB]=resampleTS_MF_const(coef, nj, B1, B2, Block, J1flag);
% Time-Scale Block BS 1d
% 
% Herwig Wendt, Lyon, 2006 - 2008

try 
    J1flag; 
catch; 
    J1flag=1; 
end;

N=nj(1);
if N==0;
    [a,b]=find(nj~=0);
    sc=b(1);
    N=nj(b)*2^(sc-1);
end
if Block>1
    % initialize blocks of indices
    bx=0:N-1;
    try % fast but memory intensive
        bx=0:N-1;
        BX=repmat(bx,Block,1);
        addit=repmat([0:Block-1]', 1,N);
        BX=BX+addit;
        BX=mod(BX,N)+1;
    catch % slow but memory save
        BX=[];bx=1:N;
        for bl=1:Block
            BX=[BX; bx];
            bx=[bx(end) bx(1:end-1)];   % create overlapping blocks
        end
    end
end


% INITIALIZE CIRCULAR BLOCK BOOTSTRAP
N_range=N;
N_BS=ceil(N/Block);
N_resample=N_BS*Block;

% % INITIALIZE COEFFICIENTS
% for j=1:length(nj);
%     DAT(j).X=NaN(1,N);
%     tempN=length(coef(j).value);
%     FACT=2^(j-1);
%     DAT(j).X(1:FACT:tempN*FACT)=coef(j).value;
% end

% INITIALIZE COEFFICIENTS - NEW (border effects !)
% does take into account shift due to border effects
% does NOT YET take into account shift due to asymmetric wavelet
% HW, Lyon, 16/09/2007
%for j=1:length(nj);
for j=max(1,J1flag):length(nj);
    DAT(j).X=NaN(1,N);
    tempN=length(coef(j).value);
    FACT=2^(j-1);
    tempL=tempN*FACT-(FACT-1);
    tempCoef=NaN(1, tempL);
    tempshift=ceil((N-tempL)/2);
    DAT(j).X(1+tempshift:FACT:tempN*FACT+tempshift)=coef(j).value;
end

%% draw indices for resamples
if Block==1
    index = fix(rand(B1,N)*N)+1 ;
    if B2>1
        index2=fix(rand(B2,B1,N)*N)+1 ;
    end
else
    for bsid=1:B1
        tempid = fix(rand(1,N_BS)*N_range)+1 ;    
        index(bsid,:)=reshape(BX(:,tempid),1,[]);
        if B2>1
            for bsid2=1:B2
                tempid2 = fix(rand(1,N_BS)*N_BS)+1 ;  
                index2(bsid2,bsid,:)=reshape(BX(:,tempid(tempid2)),1,[]);
            end
        end
    end
end

for j=max(1,J1flag):length(nj);
    for bsid=1:B1
        % get the coefficients at the indices for each resample
        tempDAT=(DAT(j).X(index(bsid,:)));
        % sort out NaNs
        tempDAT=tempDAT(~isnan(tempDAT));
        % Catch case of empty resample (mainly last scale): use ordinary
        % bootstrap sample
        if isempty(tempDAT); 
            Ntemp=length(coef(j).value);
            idtemp = fix(rand(1,Ntemp)*Ntemp)+1 ;
            tempDAT=coef(j).value(idtemp);
        else
            Ntempref=length(coef(j).value);
            Ntemp=length(tempDAT);
            if Ntemp<Ntempref
                idtemp = fix(rand(1,(Ntempref-Ntemp))*(Ntempref))+1 ;
                tempDAT= [tempDAT coef(j).value(idtemp)];
            elseif Ntemp>Ntempref
                tempDAT=tempDAT(1:Ntempref);
            end
        end         
        coefB{j}{bsid}.t=tempDAT;
        if B2>1
            for bsid2=1:B2
                tempDAT2=(DAT(j).X(squeeze(index2(bsid2,bsid,:))));
                tempDAT2=tempDAT2(~isnan(tempDAT2));
                if isempty(tempDAT2);
                    Ntemp=length(coefB{j}{bsid}.t);
                    idtemp = fix(rand(1,Ntemp)*Ntemp)+1 ;
                    tempDAT2=coefB{j}{bsid}.t(idtemp);
                else
                    Ntempref=length(coefB{j}{bsid}.t);
                    Ntemp=length(tempDAT2);
                    if Ntemp<Ntempref
                        idtemp = fix(rand(1,(Ntempref-Ntemp))*(Ntempref))+1 ;
                        tempDAT2= [tempDAT2 coefB{j}{bsid}.t(idtemp)];
                    elseif Ntemp>Ntempref
                        tempDAT2=tempDAT2(1:Ntempref);
                    end
                end
                L1=length(tempDAT); L2=length(tempDAT2);
                if L1<L2; tempDAT2=tempDAT2(1:L1); end
                if L2<L1; tempDAT2=[tempDAT2 tempDAT(L2+1:end)]; end                    
                coefB{j}{bsid}.T(bsid2,:)=tempDAT2;
            end
        end
    end
end