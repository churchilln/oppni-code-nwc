function [logstatDyad, parBSest] = MFA_BS_statlog1d(data, methodWT, paramEST, paramBS, Jflag, T_S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [logstatDyad, parBSest] = MFA_BS_statlog1d(data, methodWT, paramEST, paramBS, Jflag, T_S);
% Calculate Structure functions using DWT & LWT
%
% Herwig Wendt, Lyon, 2006 - 2008

j1=paramEST.j1;
j2=paramEST.j2;
MomNul=paramEST.MomNul;
sym=paramEST.sym;

gamint=paramEST.gamint;
EstFun=paramEST.Fun;

Fun=bin2dec(num2str(EstFun));

if Fun~=1; % if not only Cumulants
q=paramEST.q;
else
    q=NaN;
end
if rem(Fun,2)==1 % if Cumulants
Cum=paramEST.Cum;
else
    Cum=NaN;
end

if Cum==1
    Cum=2; % Order of highest Cumulant estimate - at least 2
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PHASE 0: SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check input                                 
[a,b]=size(data);
if a>b; data=data'; end
%% Check DWT, LWT
if find(methodWT==1); DWT=1; else DWT=0; end
if find(methodWT==2); LWT=1; else LWT=0; end
if (~DWT)&(~LWT)
    disp('In MFA_BS_statlog: Nothing to do (return)');
    logstatDyad=NaN;
    parBSest=NaN;    
end

%% Read regression regions for zeta, Cumulants
if length(j1)>1; j1Z=j1(1); j1C=j1(2); else j1Z=j1; j1C=j1; end
if length(j2)>1; j2Z=j2(1); j2C=j2(2); else j2Z=j2; j2C=j2; end
if Jflag; J1BS=min(j1Z, j1C); else; J1BS=1; end;

%% Bootstrap Setup %%%%%%%%%%%%%%
%% Read in Bootstrap Parameters 
Method=paramBS.Method;
if find(Method==1); NOR=1; else NOR=0; end
if find(Method==2); BAS=1; else BAS=0; end
if find(Method==3); PER=1; else PER=0; end
if find(Method==4); STU=1; else STU=0; end
if find(Method==5); BASADJ=1; else BASADJ=0; end
if find(Method==6); PERADJ=1; else PERADJ=0; end
B1=paramBS.Nresamp1;
B2=paramBS.Nresamp2;
BlockLength=paramBS.blocklength;
BlockLengthW=BlockLength;
BlockLengthL=BlockLength;
%% Check what is to do at each scale: Bootstrap Parameters for resampling
% initialise: Do primary Bootstrap 
calcMethod=[2]; TTsave=0;     
if B1>1; doB1=1; else doB1=0; end
if B2>1; doB2=1; else doB2=0; end    
 % If no bootstrap Method demanded do nothing.
if ~(NOR|BAS|PER|STU|BASADJ|PERADJ); B1=1; doB1=0; doB2=0; end;    
if ~(STU|BASADJ|PERADJ); B2=1; doB2=0; end;
if ~doB1 doB2=0; end
% ADD: Return double bootstrap estimates
if (STU|BASADJ|PERADJ)&doB2
    calcMethod=[calcMethod 4]; TTsave=1;
end
% write bootstrap parameters
parBSest=struct('Nresamp1',B1,'Nresamp2',B2, 'blocklength', BlockLength, 'Method', calcMethod);

% if EstFun==0;
%     fhandle='EstFun_Cum';
%     paramEst=struct('fhandle', fhandle, 'param',[Cum]);
% else
%     fhandle='EstFun_MFA';
%     paramEst=struct('fhandle', fhandle, 'param',[q, Cum]);
% end

fhandle='flexEstFun_MFA';
paramEst=struct('fhandle', fhandle, 'q',q, 'Cum', Cum, 'EstFun', EstFun);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHASE 1 : Coefficients Leaders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[coef, leaders, nj] = DxLx1d(data, MomNul, gamint, sym);
% [coef, leaders, nj] = IncOscDyad1d(data, MomNul, gamint); disp('Increments...');

nbvoies=min(length(nj.W),length(nj.L));
scale=2.^(1:nbvoies);

for j=1:nbvoies
    supcoef(j)=coef(j).supcoef;
    supcoefL(j)=coef(j).supcoefL;
    supcoefnoint(j)=coef(j).supcoefnoint;
    mincoef(j)=coef(j).mincoef;
    supcoefid(j)=coef(j).supcoefid;
    supcoefidL(j)=coef(j).supcoefidL;
    mincoefid(j)=coef(j).mincoefid;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHASE 2 : Bootstrap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if T_S
    %% TIME SCALE BOOTSTRAP
    %-- 2c: Coefficients - Estimates and Bootstrap
    if DWT
        fprintf('\n DWT time scale bootstrap ...\n ');
        [estimatesW, LEstW]=resample1d_T_S(coef, nj.W, paramBS, paramEst, J1BS);
        % [estimatesW{j}, LEstW{j}]=ndimBSresample(AbsdqkW, parBSest, paramEst, TTsave);
        if J1BS>1; tmpjflag=1; else tmpjflag=0; end
        for j=1:nbvoies
            if j==min(j1Z, j1C); tmpjflag=0; end
            ESTW(:,j)=estimatesW{j}.t;      % read out estimates at this scale j
            % read out Bootstrap Variance Estimates and Bootstrap Resamples at this scale
            if doB1&~tmpjflag
                ESTBW(:,:,j)=estimatesW{j}.T;
                VESTW(:,j)=estimatesW{j}.stdt;
            end
            % read out double bootstrap estimates
            if doB2&~tmpjflag
                ESTBBW(:,:,:,j)=estimatesW{j}.TT;
            end
        end
    end

    %-- 2d: Leaders- Estimates and Bootstrap
    if LWT
        fprintf('\n LWT time scale bootstrap ...\n ');
        [estimatesL, LEstL]=resample1d_T_S(leaders, nj.L, paramBS, paramEst, J1BS);
        % [estimatesW{j}, LEstW{j}]=ndimBSresample(AbsdqkW, parBSest, paramEst, TTsave);
        if J1BS>1; tmpjflag=1; else tmpjflag=0; end
        for j=1:nbvoies
            if j==min(j1Z, j1C); tmpjflag=0; end
            ESTL(:,j)=estimatesL{j}.t;      % read out estimates at this scale j
            % read out Bootstrap Variance Estimates and Bootstrap Resamples at this scale
            if doB1&~tmpjflag
                ESTBL(:,:,j)=estimatesL{j}.T;
                VESTL(:,j)=estimatesL{j}.stdt;
            end
            % read out double bootstrap estimates
            if doB2&~tmpjflag
                ESTBBL(:,:,:,j)=estimatesL{j}.TT;
            end
        end
    end
else
    %% ORDINARY BOOTSTRAP
    BLXW = BlockLengthW*ones(1, nbvoies);
    BLXL = BlockLengthL*ones(1, nbvoies);
    for j=1:nbvoies
       %-- 2a: Check block size
        if DWT
            while (nj.W(j) < 2*BlockLengthW) & (BlockLengthW~=1) % check if block size is appropriate
                BlockLengthW = max(fix(BlockLengthW/2),1);  BLXW(j:nbvoies) = BlockLengthW*ones(size(j:1:nbvoies)) ;
            end
        end
        if LWT
            while (nj.L(j) < 2*BlockLengthL)  & (BlockLengthL~=1)
                BlockLengthL = max(fix(BlockLengthL/2),1);  BLXL(j:nbvoies) = BlockLengthL*ones(size(j:1:nbvoies)) ;
            end
        end
        if DWT&LWT
            parBSest.blocklength=min(BlockLengthW,BlockLengthL);
        elseif DWT
            parBSest.blocklength=BlockLengthW;
        elseif LWT
            parBSest.blocklength=BlockLengthL;
        end

        %-- 2b: Check if Bootstrap
        if j==min(j1Z, j1C); Jflag=0;  parBSest.Nresamp1=B1; parBSest.Nresamp2=B2; end; % Reset Jflag to 0 for j>=j1: do bootstrap
        if Jflag    % don't do bootstrap at this scale
            parBSest.Nresamp1=1;
            parBSest.Nresamp2=1;
        end

        %-- 2c: Coefficients - Estimates and Bootstrap
        if DWT
            [estimatesW{j}, LEstW{j}]=resample1d(coef(j).value, parBSest, paramEst, TTsave);
            ESTW(:,j)=estimatesW{j}.t;      % read out estimates at this scale j
            % read out Bootstrap Variance Estimates and Bootstrap Resamples at this scale
            if doB1&~Jflag
                ESTBW(:,:,j)=estimatesW{j}.T;
                VESTW(:,j)=estimatesW{j}.stdt;
            end
            % read out double bootstrap estimates
            if doB2&~Jflag
                ESTBBW(:,:,:,j)=estimatesW{j}.TT;
            end
        end

        %-- 2d: Leaders- Estimates and Bootstrap
        if LWT
            [estimatesL{j}, LEstL{j}]=resample1d(leaders(j).value, parBSest, paramEst, TTsave);
            ESTL(:,j) =estimatesL{j}.t; % read out estimates at this scale j
            % read out Bootstrap Variance Estimates and Bootstrap Resamples at this scale
            if doB1&~Jflag
                ESTBL(:,:,j)=estimatesL{j}.T;
                VESTL(:,j)=estimatesL{j}.stdt;
            end
            % read out double bootstrap estimates
            if doB2&~Jflag
                ESTBBL(:,:,:,j)=estimatesL{j}.TT;
            end
        end        
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHASE 3 : Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j2=min(j2, nbvoies);
% output Coefficients
if DWT
logstatDyad.DWT.gamint=gamint;
logstatDyad.DWT.est=ESTW;          
logstatDyad.DWT.supcoef=supcoef;
logstatDyad.DWT.supcoefL=supcoefL;
logstatDyad.DWT.supcoefnoint=supcoefnoint;
logstatDyad.DWT.mincoef=mincoef;
logstatDyad.DWT.supcoefid=supcoefid;
logstatDyad.DWT.supcoefidL=supcoefidL;
logstatDyad.DWT.mincoefid=mincoefid;
logstatDyad.DWT.nj=nj.W;             
logstatDyad.DWT.Lest=LEstW;       
logstatDyad.DWT.EstFun=EstFun;
if doB1
    logstatDyad.DWT.estB=ESTBW;     
    logstatDyad.DWT.Vest=VESTW;     
end
if doB2
    logstatDyad.DWT.estBB=ESTBBW;   
end
logstatDyad.DWT.scale=scale;
logstatDyad.DWT.j1=j1;
logstatDyad.DWT.j2=j2;
end
% output Leaders
if LWT
logstatDyad.LWT.gamint=gamint;
logstatDyad.LWT.est=ESTL;          
logstatDyad.LWT.supcoef=supcoef;
logstatDyad.LWT.supcoefL=supcoefL;
logstatDyad.LWT.supcoefnoint=supcoefnoint;
logstatDyad.LWT.mincoef=mincoef;
logstatDyad.LWT.supcoefid=supcoefid;
logstatDyad.DWT.supcoefidL=supcoefidL;
logstatDyad.LWT.mincoefid=mincoefid;
logstatDyad.LWT.nj=nj.L;             
logstatDyad.LWT.Lest=LEstL;    
logstatDyad.LWT.EstFun=EstFun;
if doB1
    logstatDyad.LWT.estB=ESTBL;     
    logstatDyad.LWT.Vest=VESTL;     
end
if doB2
    logstatDyad.LWT.estBB=ESTBBL;   
end
logstatDyad.LWT.scale=scale;
logstatDyad.LWT.j1=j1;
logstatDyad.LWT.j2=j2;
end
parBSest=struct('Nresamp1',B1,'Nresamp2',B2, 'blocklength', BlockLength, 'Method', Method, 'T_S', T_S);
