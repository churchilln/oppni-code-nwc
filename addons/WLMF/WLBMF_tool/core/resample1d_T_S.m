function [estimates, LEst]=resample1d_T_S(coef, nj, paramBS, paramEst, J1flag);
% function [estimates, LEst]=resample1d_T_S(coef, nj, paramBS, paramEst, J1flag);
% Time-Scale Block BS 1d
% 
% Herwig Wendt, Lyon, 2006 - 2008

try 
    J1flag; 
catch; 
    J1flag=1; 
end;

InputError1='paramBS must be a structure with elements:\n   paramBS.Nresamp1 \n   paramBS.Nresamp2\n   paramBS.blocklength\n   paramBS.Method\n   paramBS.verb (optional)\n';
try  NB1=isreal(paramBS.Nresamp1);      if ~NB1; paramBS.Nresamp1=99; end;  catch error(InputError1, NB1); end
try  NB2=isreal(paramBS.Nresamp2);      if ~NB2; paramBS.Nresamp2=1;  end; catch error(InputError1, NB2); end
try  NB3=isreal(paramBS.blocklength);   if ~NB3; paramBS.blocklength=floor(nj(1)/32); end;  catch error(InputError1, NB3); end
try  NB4=isreal(paramBS.Method);        if ~NB4; paramBS.Method=[1:6];end;  catch error(InputError1, NB4); end
try CHR=ischar(paramEst.fhandle); catch error('The structure paramEst must contain a field  paramEst.fhandle with a valid function handle'); end
if ~CHR; error('The function handle paramEst.fhandle is not valid'); end

fhandle=paramEst.fhandle;   % Estimator name
FHandle=str2func(fhandle);       % Function handle for Estimator
NOUT=nargout(FHandle);           % Number of Output Arguments of Estimator
% Check if Estimator needs extra parameters
% fparam=1;
% if length(fieldnames(paramEst))==1
%     fparam=0;
% elseif isempty(paramEst.param)
%     fparam=0;
% end

EstFun=paramEst.EstFun;
Fun=bin2dec(num2str(EstFun));

if (Fun==1)||(Fun==4); NOUT=1;
elseif (Fun==2)||(Fun==5); NOUT=2;
elseif (Fun==3)||(Fun==6); NOUT=3;
else; NOUT=4; end

% Check which Bootstrap statistics are to be calculated
Method=paramBS.Method; CIdisp=[]; Hstat=[];
if find(Method==1); NOR=1; else NOR=0; end
if find(Method==2); BAS=1; else BAS=0; end
if find(Method==3); PER=1; else PER=0; end
if find(Method==4); STU=1; else STU=0; end
if find(Method==5); BASADJ=1; else BASADJ=0; end
if find(Method==6); PERADJ=1; else PERADJ=0; end

if STU|PERADJ|BASADJ; doB2=1; else; doB2=0; end

B1=paramBS.Nresamp1;        % number of primary bootstrap resamples
B2=paramBS.Nresamp2;        % number of bootstrap resamples for variance estimates (used for Normal, Studentised, Variance Stabilising Transformation)
Block=paramBS.blocklength;  % Block length for moving blocks resampling

N=nj(1);
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
    if doB2;
        indexB2 = fix(rand(B2,B1,N)*N)+1 ;
    end
else
    for bsid=1:B1
        tempid = fix(rand(1,N_BS)*N_range)+1 ;    
        index(bsid,:)=reshape(BX(:,tempid),1,[]);
        if doB2;
            for bs2id=1:B2
                tempidB2 = fix(rand(1,N_BS)*N_BS)+1 ;    
                indexB2(bs2id,bsid,:)=reshape(BX(:,tempid(tempidB2)),1,[]); 
            end
        end
    end
end

for j=1:length(nj);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Structure Functions
    t=cell(1,NOUT);
%     if ~fparam
%         %[t_out{:}]=feval(fhandle,X);
%         [t{:}]=FHandle(coef(j).value);
%     else
        %[t{:}]=FHandle(coef(j).value,paramEst.param);
        [t{:}]=FHandle(coef(j).value,paramEst);
%     end
    estimates{j}.t=cell2mat(t);
    % get length of estimates
    for ii=1:NOUT
        LEst{j}(ii)=length(t{ii});
    end
    
    if j>=J1flag
        T=cell(B1,NOUT);
        if doB2; TT=cell(B2,B1,NOUT); end;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            % Calculate BS Structure Functions            
%             if ~fparam
%                 [T{bsid,:}]=FHandle(tempDAT);
%             else
                %[T{bsid,:}]=FHandle(tempDAT, paramEst.param);
                [T{bsid,:}]=FHandle(tempDAT, paramEst);
%             end
            
            if doB2;
                for bs2id=1:B2
                    tmpTT=cell(1,NOUT);
                    % get the coefficients at the indices for each resample
                    tempDAT2=(DAT(j).X(squeeze(indexB2(bs2id,bsid,:))));
                    % sort out NaNs
                    tempDAT2=tempDAT2(~isnan(tempDAT2));
                    % Catch case of empty resample (mainly last scale): use ordinary
                    % bootstrap sample                    
                    if isempty(tempDAT2);
                        Ntemp=length(tempDAT);
                        idtemp = fix(rand(1,Ntemp)*Ntemp)+1 ;
                        tempDAT2=tempDAT(idtemp);
                    else
                        Ntempref=length(tempDAT);
                        Ntemp=length(tempDAT2);
                        if Ntemp<Ntempref
                            idtemp = fix(rand(1,(Ntempref-Ntemp))*(Ntempref))+1 ;
                            tempDAT2= [tempDAT2 tempDAT(idtemp)];
                        elseif Ntemp>Ntempref
                            tempDAT2=tempDAT2(1:Ntempref);
                        end
                    end
                    % Calculate BS Structure Functions
%                     if ~fparam
%                         [tmpTT{:}]=FHandle(tempDAT2);
%                     else
                        %[tmpTT{:}]=FHandle(tempDAT2, paramEst.param);
                        [tmpTT{:}]=FHandle(tempDAT2, paramEst);
%                     end
                    tempTT(bs2id,bsid,:)=cell2mat(tmpTT);
               end
            end
        end
        estimates{j}.T=cell2mat(T);
        estimates{j}.stdt=std(estimates{j}.T);

        if doB2; 
            estimates{j}.stdT=squeeze(std(tempTT));
            estimates{j}.TT=tempTT; 
        end
    else
        estimates{j}.T=NaN(B1, length(estimates{j}.t));
        estimates{j}.stdt=NaN(1, length(estimates{j}.t));
        if doB2; 
            estimates{j}.stdT=NaN(B1, length(estimates{j}.t));
            estimates{j}.TT=NaN(B2, B1, length(estimates{j}.t)); 
        end
    end
end

    
