function [incr, osc, nj] = IncOscDyad1d(data, Nwt, gamint)

if nargin<3; gamint=0; end;

hh0=[-1 1]/2; % first order increment (L1 Norm)
hh=1;
for nn=1:Nwt;
    hh=conv(hh,hh0); % Nwt-th order increment (L1 Norm)
end
nl=length(hh);

%--- Predict the max # of octaves available given Nwt
n = length(data) ;                   % data length
nbvoies= fix( log2(length(data)) );
nbvoies = min( fix(log2(n/(2*Nwt+1)))  , nbvoies);  %   safer, casadestime having problems

njtemp = n;    %  this is the the appro, changes at each scale
appro=data;
for j=1:nbvoies         % Loop Scales
    %-- Phase 1a: get the increment coefficients at this scale
    njtemp = length(appro) ;
    conv_hh = conv(appro,hh) ;
    appro = appro(1:2:njtemp) ; %-- prepare for next coarser scale
    decime = conv_hh(nl:2:end-nl+1) ;   % details at scale j. decime always becomes empty near the end,
    nj.W(j) = length(decime);            % number of coefficients
%     incr(j).value_noabs=decime;
%     decime = abs(decime);
%     incr(j).value=decime;
    
    AbsdqkW = abs(decime);   %%%%%%% passage Norme L1
    %% ----> ADDED 16/04/2007
    % hmin before integration
    incr(j).supcoefnoint=max(AbsdqkW);
    % fractional integration
    AbsdqkW = AbsdqkW*2^(gamint*j);
    [incr(j).supcoef, incr(j).supcoefid]=max(AbsdqkW);
    incr(j).value=AbsdqkW;%(find(AbsdqkW>=eps));
    incr(j).value_noabs=decime*2^(gamint*j);
    incr(j).gamma=gamint;
    incr(j).sign=sign(decime);
    
    decime = abs(decime)*2^(gamint*j);
    %-- Phase 1b: calculate the oscillations at this scale
    if j==1
        osc(1).value = decime;
        osc(1).value3lam = max([decime(1:end-2);decime(2:end-1);decime(3:end)]);
        nj.L(j) = length(osc(1).value3lam);
    else
        length_detail = 2*length(decime);
        osc(j).value = max([decime; osc(j-1).value(1:2:length_detail); osc(j-1).value(2:2:length_detail)]);
        osc(j).value3lam = max([osc(j).value(1:end-2);osc(j).value(2:end-1);osc(j).value(3:end)]);
        nj.L(j) = length(osc(j).value3lam);
        if nj.L(j) == 0                      % forget this j and exit if no coefficients left!
            fprintf('In wtspecq_statlog, oops!  no details left at scale %d! \n',j)
        end
    end
    
    osc(j).gamma=gamint;    
    [osc(j).mincoef, osc(j).mincoefid]=min(osc(j).value);
    [osc(j).supcoefL, osc(j).supcoefidL]=max(osc(j).value);
    osc(j).supcoefnoint=incr(j).supcoefnoint; osc(j).supcoef=incr(j).supcoef; osc(j).supcoefid=incr(j).supcoefid;
    incr(j).mincoef=osc(j).mincoef; incr(j).mincoefid=osc(j).mincoefid;
    incr(j).supcoefL=osc(j).supcoefL; incr(j).supcoefidL=osc(j).supcoefidL;
end