function [est_m, logstat_m, est_all] = MF_const_BS(coef, nj, M, paramBS, paramEst);
% function [est_m, logstat_m, est_all] = MF_const_BS(coef, nj, M, paramBS, paramEst);
% calculate bootstrap resamples for blocks of time series
%
% Herwig Wendt, Lyon, 2006 - 2008

loge=log2(exp(1));

% Get Parameters
B1=paramBS.Nresamp1;
B2=paramBS.Nresamp2;
Block=paramBS.blocklength;
Method=paramBS.Method;
T_S=paramBS.T_S;
j1=paramEst.j1;
j2=paramEst.j2;
wtype=paramEst.wtype;
Jflag=paramEst.Jflag;
Cum=paramEst.Cum;
q=paramEst.q;
EstFun=paramEst.EstFun;

Fun=0;
    if EstFun>=100; Fun=Fun+4; EstFun=EstFun-100; end
    if EstFun>=10; Fun=Fun+2; EstFun=EstFun-10; end
    Fun=Fun+EstFun;
    
poscount=1;
if Fun>=4;
    ZQ=1;
    zqpos=1; poscount=poscount+1; 
else
    ZQ=0;
end
if (Fun~=1) && (Fun~=4) && (Fun~=5);
    Dqpos=poscount;  poscount=poscount+1;
    hqpos=poscount;  poscount=poscount+1;
    DH=1;
else
    DH=0;
end
if rem(Fun,2)==1  
    Cppos=poscount;
    CP=1;
else
    CP=0;
end
NOUT=ZQ+2*DH+CP;

if (NOUT==1)&((Cum==1)|length(q)==1); Lest1=1; else; Lest1=0; end

FHandle=str2func('flexEstFun_MFA');

if T_S; j2max=find(nj>=8*M); else; j2max=find(nj>=4*M); end
j2max=j2max(end); j2=min(j2, j2max);
nD=floor(nj/M);
%posC=floor(nj/2);
if Jflag; Jstart=j1; else; Jstart=1; end

%% whole series: Estimates
if nargout>2
    j2all=j2+floor(log2(M));
    % --- get structure functions
    for j=Jstart:j2all
        Dxtmp=coef(j).value;
        t=cell(1,NOUT);
        [t{:}]=FHandle(Dxtmp,paramEst);
        for ii=1:NOUT
            LEstall(ii)=length(t{ii});
        end
        Sqjall(:,j)=cell2mat(t);
    end
    LEWall=cumsum([0 LEstall]);
    if CP; Sqjall(LEWall(Cppos)+1:LEWall(Cppos+1),:)=Sqjall(LEWall(Cppos)+1:LEWall(Cppos+1),:)*loge; end
    % --- get regressions
%    if ~Lest1
        [est_all.t,V,Q,aest]=MFA_BS_regrmat(Sqjall,0, nj, wtype, j1, j2all);
        if DH;
            est_all.t(LEWall(Dqpos)+1:LEWall(Dqpos+1))=est_all.t(LEWall(Dqpos)+1:LEWall(Dqpos+1))+1;
        end
%    else
%        [est_all.t,Vzeta,Q,aest]=regrespond_det2(2,3,Sqjall ,ones(size(Sqjall)),nj, wtype,  j1,j2,0);
%    end    
    est_all.nj=nj;
    est_all.j1=j1;
    est_all.j2=j2all;
end
%% subseries: estimates and bootstrap estimates
clear coefB coef_m coefB_m coefBsig_m;
if B2<=1;    B2sig=B1; else;    B2sig=B2; end

% --- LEVEL 2: Bootstrap resamples from whole series
if T_S
    [coefB]=resampleTS_MF_const(coef, nj, B1, B2, Block, Jstart);
else
    for j=Jstart:j2
        [coefB{j}]=resample_MF_const(coef(j).value,Block,B1,B2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- BLOCKS
for m=1:M
    clear coef_m coefBsig_m coefB_m
    % --- LEVEL 2 : Cut whole series BS resamples of coefficients

    for j=Jstart:j2
        for b=1:B1
            coefB_m{j}{b}.t=coefB{j}{b}.t((m-1)*nD(j)+1:min(m*nD(j),length(coefB{j}{b}.t)));
            if B2>1;
                coefB_m{j}{b}.T=coefB{j}{b}.T(:,(m-1)*nD(j)+1:min(m*nD(j),length(coefB{j}{b}.t)));
            end
        end
    end
   
    % --- LEVEL 1 : create blocks of coefficients
    %               bootstrap resample for variance estimate
    for j=1:j2
        coef_m(j).value=coef(j).value((m-1)*nD(j)+1:min(m*nD(j),nj(j)));
        nj_m(m,j)=length(coef_m(j).value);
        if (~T_S)&(j>=Jstart)
            [coefBsig_m{j}]=resample_MF_const(coef_m(j).value,Block,B2sig, 1);
        end
    end
    if T_S
        [coefBsig_m]=resampleTS_MF_const(coef_m, nj_m(m,:), B2sig, 1, Block, Jstart);
    end    
    
    % ------- Get Structure Function
    for j=Jstart:j2
        if ~((j<j1)&Jflag)
            % ---  LEVEL 1
            % Estimates level 1
            Dx{m}{j}=coef_m(j).value;
            t=cell(1,NOUT);
            [t{:}]=FHandle(Dx{m}{j},paramEst);
            for ii=1:NOUT
                LEst(ii)=length(t{ii});
            end                       
            Sqj{m}(:,j)=cell2mat(t);
            % Bootstrap estimates level 1
            if B2sig>1
                t=cell(1,NOUT);
                for b=1:B2sig
                    DxBsig{m}{j}=coefBsig_m{j}{b}.t;                    
                    [t{:}]=FHandle(DxBsig{m}{j},paramEst);                    
                    SqjBsig{m}(b,:,j)=cell2mat(t);
                end
            end
            % --- LEVEL 2
            % Bootstrap Estimates level 2
            % Double Bootstrap Estimates level 2
            for b=1:B1
                DxB{m}{j}=coefB_m{j}{b}.t;
                [t{:}]=FHandle(DxB{m}{j},paramEst);
                SqjB{m}(b,:,j)=cell2mat(t);
                if B2>1;                    
                    DxBB{j}=coefB_m{j}{b}.T;
                    for b2=1:B2
                        [t{:}]=FHandle(DxBB{j}(b2,:),paramEst);
                        SqjBB{m}(b2, b,:,j)=cell2mat(t);
                    end
                end
            end  % for b=1:B1
        end % if Jflag
    end % for j= ...
    LEW=cumsum([0 LEst]);
    if CP
        Sqj{m}(LEW(Cppos)+1:LEW(Cppos+1),:)=Sqj{m}(LEW(Cppos)+1:LEW(Cppos+1),:)*loge;
        SqjB{m}(:,LEW(Cppos)+1:LEW(Cppos+1),:)=SqjB{m}(:,LEW(Cppos)+1:LEW(Cppos+1),:)*loge;
        if B2>1;
            SqjBB{m}(:,:,LEW(Cppos)+1:LEW(Cppos+1),:)=SqjBB{m}(:,:,LEW(Cppos)+1:LEW(Cppos+1),:)*loge;
        end
        if B2sig>1
            SqjBsig{m}(:,LEW(Cppos)+1:LEW(Cppos+1),:)=SqjBsig{m}(:,LEW(Cppos)+1:LEW(Cppos+1),:)*loge;
        end
    end

    % ------- Get Regressions
    % get regression range
    est_m{m}.j1=j1; [tmp, J2]=size(Sqj{m}); est_m{m}.j2=min(J2,j2); j2=est_m{m}.j2;

    %calculate theoretical variance of zeta(2)
    jj = j1:1:j2 ;J = length(jj); loge=log2(exp(1)); precis_zeta= 10^-5;
    warning off;
    varj = 2./nj_m(m,:)*loge^2;
    warning on;
    varjj = varj(jj) ;yjj=Sqj{m}(jj);
    S0 = sum(1./varjj) ;
    S1 = sum(jj./varjj) ;
    S2 = sum(jj.^2./varjj) ;
    wjj = (S0 * jj - S1) ./ varjj / (S0*S2-S1*S1) ;
    vjj = (S2 - S1 * jj) ./ varjj / (S0*S2-S1*S1) ;
    alphaest  = sum(wjj .* yjj ) ;
    Valpha = sum(varjj.*wjj.*wjj) ;
    est_m{m}.Vtheo=Valpha;

    % Block Estimates
    [est_m{m}.t,V,Q,aest]=MFA_BS_regrmat(Sqj{m},0, nj_m(m,:), wtype, j1, j2);
    % Bootstrap Block Estimates
    SqjB{m}=squeeze(SqjB{m});
    
    [est_m{m}.T,V,Q,aest]=MFA_BS_regrmat(SqjB{m},0,nj_m(m,:), wtype, j1, j2);
    if DH;
        est_m{m}.t(LEW(Dqpos)+1:LEW(Dqpos+1))=est_m{m}.t(LEW(Dqpos)+1:LEW(Dqpos+1))+1;
        est_m{m}.T(LEW(Dqpos)+1:LEW(Dqpos+1),:)=est_m{m}.T(LEW(Dqpos)+1:LEW(Dqpos+1),:)+1;
    end
    if B2>1
        % Bootstrap and double bootstrap STD estimation
        SqjBB{m}=squeeze(SqjBB{m});
        SqjBsig{m}=squeeze(SqjBsig{m});
        [estsig_m{m}.T,V,Q,aest]=MFA_BS_regrmat(SqjBsig{m},0,nj_m(m,:), wtype, j1, j2);        
        [TTtmp,V,Q,aest]=MFA_BS_regrmat(SqjBB{m},0,nj_m(m,:), wtype, j1, j2);
        if DH;
            estsig_m{m}.T(LEW(Dqpos)+1:LEW(Dqpos+1),:)=estsig_m{m}.T(LEW(Dqpos)+1:LEW(Dqpos+1),:)+1;
            TTtmp(:,LEW(Dqpos)+1:LEW(Dqpos+1),:)=TTtmp(:,LEW(Dqpos)+1:LEW(Dqpos+1),:)+1;
        end
        est_m{m}.stdT=std(estsig_m{m}.T, 0, 2)';
        est_m{m}.stdTT=std(TTtmp, 0, 3-Lest1)';
    elseif B2sig>0
        % Bootstrap STD estimation, used as well as double bootstrap STD
        SqjBsig{m}=squeeze(SqjBsig{m});
        [estsig_m{m}.T,V,Q,aest]=MFA_BS_regrmat(SqjBsig{m},0,nj_m(m,:), wtype, j1, j2);
        if DH;
            estsig_m{m}.T(LEW(Dqpos)+1:LEW(Dqpos+1),:)=estsig_m{m}.T(LEW(Dqpos)+1:LEW(Dqpos+1),:)+1;
        end
        est_m{m}.stdT=std(estsig_m{m}.T, 0, 2)';
        est_m{m}.stdTT=repmat(std(est_m{m}.T, 0, 2),1,B1);
    end
    
    % get data length
    est_m{m}.n=length(coef_m(1).value)*2;
    est_m{m}.LEst=LEW;
    est_m{m}.nj=nj_m(m,:);
    logstat_m{m}.Sqj=Sqj{m};
    logstat_m{m}.SqjBcut=SqjBsig{m};
    logstat_m{m}.SqjBall=SqjB{m};
end % loop on m
