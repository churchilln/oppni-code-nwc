function [estDyad]=MFA_BS_regrest(logstat, paramBS, paramEST, wtype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [estDyad]=MFA_BS_regrest(logstat, paramBS, paramEST, wtype)
%   calculate final estimates from structure functions
% 
% Herwig Wendt, Lyon, 2006 - 2008

EstFun=paramEST.Fun;
Fun=bin2dec(num2str(EstFun));

loge=log2(exp(1));
try logstat.regrmult; regcor=logstat.regrmult; catch; regcor=1; end
try logstat.j1; j1=logstat.j1; catch;  end
try logstat.j2; j2=logstat.j2; catch;  end
%% Read regresion regions for zeta, Cumulants
if length(j1)>1; j1Z=j1(1); j1C=j1(2); else j1Z=j1; j1C=j1; end
if length(j2)>1; j2Z=j2(1); j2C=j2(2); else j2Z=j2; j2C=j2; end
LEW=[0 cumsum(logstat.Lest{1})];       % Length of estimates
%Cum=LEW(end)-LEW(end-1);

poscount=1;
if Fun>=4;
    zqpos=1; poscount=poscount+1;    
end
if (Fun~=1) && (Fun~=4) && (Fun~=5);
    Dqpos=poscount;  poscount=poscount+1;
    hqpos=poscount;  poscount=poscount+1;
    DH=1;
else
    DH=0;
end
if rem(Fun,2)==1
    Cum=paramEST.Cum;    
    Cppos=poscount;
    CP=1;
else
    CP=0;
end

try logstat.imagedata; dimcor=logstat.imagedata+1; catch; dimcor=1; end
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
if B1>1; doB1=1; else doB1=0; end
if B2>1; doB2=1; else doB2=0; end 
if ~(NOR|BAS|PER|STU|BASADJ|PERADJ); B1=1; doB1=0; doB2=0; end;    
if ~(STU|BASADJ|PERADJ); B2=1; doB2=0; end;
if ~doB1; doB2=0; end

%%%%% ADDED 08/12/06 %%%%%
%NEW=1;
%if NEW
regcor=1;
%end;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in Estimates
% Coefficients
ESTW=logstat.est;          
njW=logstat.nj;         
LEstW=logstat.Lest;        
if doB1
    ESTBW=logstat.estB;     
    VESTW=logstat.Vest;     
end
if doB2
    ESTBBW=logstat.estBB;   
end

%% --- > UNIFORM REGULARITY H_min
[h_min,Vzeta,QW,h_min_aest]=MFA_BS_regrmat(log2(logstat.supcoef),ones(size(logstat.supcoef)),njW, wtype, j1C,j2C);
[h_minL,Vzeta,QW,h_minL_aest]=MFA_BS_regrmat(log2(logstat.supcoefL),ones(size(logstat.supcoefL)),njW, wtype, j1C,j2C);
[h_minnoint,Vzeta,QW,h_minnoint_aest]=MFA_BS_regrmat(log2(logstat.supcoefnoint),ones(size(logstat.supcoefnoint)),njW, wtype, j1C,j2C);
%% < --- 
%% --- > UNIFORM REGULARITY H_max
[h_max,Vzeta,QW,h_max_aest]=MFA_BS_regrmat(log2(logstat.mincoef),ones(size(logstat.mincoef)),njW, wtype, j1C,j2C);
%% < --- 


% catch case of no bootstrap
if ~(B1>1); 
    VESTW=ones(size(ESTW)); %VESTL=VESTW;
end
%-- 3a: ESTIMATES
% Normalize Cumulants(j) for regression
% if Nest>1
%     ESTW((LEW(4)+1):LEW(5),:)=ESTW((LEW(4)+1):LEW(5),:)*loge;    
% else
%     ESTW=ESTW*loge;
% end
if CP
    ESTW((LEW(Cppos)+1):LEW(Cppos+1),:)=ESTW((LEW(Cppos)+1):LEW(Cppos+1),:)*loge; 
end
% Regression - cum
%if NEW
    [REstW,Vzeta,QW,aestW]=MFA_BS_regrmat(ESTW,VESTW.^2,njW, wtype, j1C,j2C);
%else
%    [REstW,Vzeta,QW,aestW]=regrmat(ESTW,VESTW.^2,njW, wtype,  j1C,j2C,0, 'W:');
%end
REstW=REstW*regcor;
% if Nest>1
%     REstW(LEW(2)+1:LEW(3))=REstW(LEW(2)+1:LEW(3))+dimcor;
% end
if DH
    REstW(LEW(Dqpos)+1:LEW(Dqpos+1))=REstW(LEW(Dqpos)+1:LEW(Dqpos+1))+dimcor;
end
if regcor~=1
    aestW= aestW-REstW*(log2(logstat.scale(1))-1/regcor);
end
%-- 3b: BOOTSTRAP ESTIMATES
if doB1
    % reshape variance for n-dim regression 
    VVESTW=squeeze(repmat(VESTW,[1 1 B1]));
    VVESTW=shiftdim(VVESTW,length(size(VVESTW))-1);
    % Normalize Cumulants(j) for regression
%     if Nest>1
%         ESTBW(:,(LEW(4)+1):LEW(5),:)=ESTBW(:,(LEW(4)+1):LEW(5),:)*loge;
%     else
%         ESTBW=ESTBW*loge;
%     end
    if CP;
        ESTBW(:,(LEW(Cppos)+1):LEW(Cppos+1),:)=ESTBW(:,(LEW(Cppos)+1):LEW(Cppos+1),:)*loge;
    end
    % Regression - cum
%    if NEW
        [REstBW,Vzeta,Q,aest]=MFA_BS_regrmat(ESTBW,VVESTW.^2,njW, wtype, j1C,j2C);
%    else
%        [REstBW ,Vzeta,Q,aest]=regrmat(ESTBW,VVESTW.^2,njW, wtype,  j1C,j2C,0);
%    end
    % Put in correct shape   
    REstBW=REstBW';
    REstBW=REstBW*regcor;
%     if Nest>1
%         REstBW(:,LEW(2)+1:LEW(3))=REstBW(:,LEW(2)+1:LEW(3))+dimcor;
%     end
    if DH
        REstBW(:,LEW(Dqpos)+1:LEW(Dqpos+1))=REstBW(:,LEW(Dqpos)+1:LEW(Dqpos+1))+dimcor;
    end
	% Bootstrap Standard Deviation
        sigmaREstBW = std(REstBW);
end

%-- 3c: DOUBLE BOOTSTRAP ESTIMATES
if doB2
    % reshape variance for n-dim regression 
    VVVESTW=squeeze(repmat(VVESTW,[1 1 1 B2]));
    VVVESTW=shiftdim(VVVESTW,length(size(VVVESTW))-1);
    % Normalize Cumulants(j) for regression
%     if Nest>1
%         ESTBBW(:,:,(LEW(4)+1):LEW(5),:)=ESTBBW(:,:,(LEW(4)+1):LEW(5),:)*loge;
%     else
%         ESTBBW=ESTBBW*loge;
%     end
    if CP
        ESTBBW(:,:,(LEW(Cppos)+1):LEW(Cppos+1),:)=ESTBBW(:,:,(LEW(Cppos)+1):LEW(Cppos+1),:)*loge;
    end
    % Regression - cum
%    if NEW
        [REstBBW,Vzeta,Q,aest]=MFA_BS_regrmat(ESTBBW,VVVESTW.^2,njW, wtype,  j1C,j2C);
%    else
%        [REstBBW,Vzeta,Q,aest]=regrmat(ESTBBW,VVVESTW.^2,njW, wtype,  j1C,j2C,0);
%    end
    % Put in correct shape
    REstBBW=permute(REstBBW, [3 1 2]);
    REstBBW=REstBBW*regcor;
%     if Nest>1
%         REstBBW(:,:,LEW(2)+1:LEW(3))=REstBBW(:,:,LEW(2)+1:LEW(3))+dimcor;
%     end
    if DH
        REstBBW(:,:,LEW(Dqpos)+1:LEW(Dqpos+1))=REstBBW(:,:,LEW(Dqpos)+1:LEW(Dqpos+1))+dimcor;
    end
    % Double Bootstrap Standard Deviation
        sigmaREstBBW = squeeze(std(REstBBW));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Phase 4: OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~(NOR|BAS|PER|STU|BASADJ|PERADJ)|~doB1; 
    estW=struct('t', REstW, 'LEst', LEstW{end});
end

if (NOR|BAS|PER)&doB1
    estW=struct('t', REstW, 'LEst', LEstW{end}, 'stdt', sigmaREstBW, 'T', REstBW);
end

if (STU|BASADJ|PERADJ)&doB2
    estW=struct('t', REstW, 'LEst', LEstW{end}, 'stdt', sigmaREstBW, 'T', REstBW, 'stdT', sigmaREstBBW);
    if BASADJ|PERADJ
        estW=struct('t', REstW, 'LEst', LEstW{end}, 'stdt', sigmaREstBW, 'T', REstBW, 'stdT', sigmaREstBBW, 'TT', REstBBW);
    end
end
estW.aest=aestW;
estW.Q=QW;
estW.h_min=h_min;
estW.h_min_aest=h_min_aest;
estW.h_minL=h_minL;
estW.h_minL_aest=h_minL_aest;
estW.h_minnoint=h_minnoint;
estW.h_minnoint_aest=h_minnoint_aest;
estW.h_max=h_max;
estW.h_max_aest=h_max_aest;
estDyad=estW;
