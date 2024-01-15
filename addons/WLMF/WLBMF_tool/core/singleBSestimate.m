function [estsingle] = singleBSest(EST, estid, paramBS);
% function [estsingle] = singleBSest(EST, estid, paramBS);
% read out single (bootstrap) estimates from n dim estimate obtained by
% MFA_regrest_BS: for BSlimit and BShtest
%
% Herwig Wendt, Lyon, 2006 - 2008

%% check if estimate is available
if estid>length(EST.t)
    error('Estimate estid not available');
end

%% check which BS estimates are available and necessary
Method = paramBS.Method;
    if find(Method==1); NOR=1; else NOR=0; end
    if find(Method==2); BAS=1; else BAS=0; end
    if find(Method==3); PER=1; else PER=0; end
    if find(Method==4); STU=1; else STU=0; end
    if find(Method==5); BASADJ=1; else BASADJ=0; end
    if find(Method==6); PERADJ=1; else PERADJ=0; end

%% write appropriate structure
if NOR|BAS|PER
    estsingle=struct('t', EST.t(estid), 'stdt', EST.stdt(estid), 'T', EST.T(:,estid)');
end
if (STU|BASADJ|PERADJ)
    estsingle=struct('t', EST.t(estid), 'stdt', EST.stdt(estid), 'T', EST.T(:,estid)', 'stdT', EST.stdT(:,estid)');
    if BASADJ|PERADJ
        estsingle=struct('t', EST.t(estid), 'stdt', EST.stdt(estid), 'T', EST.T(:,estid)', 'stdT', EST.stdT(:,estid)', 'TT', EST.TT(:,:,estid)');
    end
end