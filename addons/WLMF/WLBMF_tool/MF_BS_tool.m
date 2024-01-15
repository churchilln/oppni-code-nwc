function [varargout] = MF_BS_tool(data, paramEST, methodWT, paramBS, paramTest, verbose, FigNum);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [varargout] = MFA_BS_light(data, paramEST, methodWT, paramBS, paramTest, verbose, FigNum);
% - Calculate cp, (zeta(q), D(h))
% - Bootstrap CI
% - Bootstrap tests
% Analysis Methods: DWT, LWT
%
% SEE demo_MF_BS_tool for usage, outputs and parameters
%
% Herwig Wendt, Lyon, 2006 - 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use the toolbox, please quote:
%
%     @ARTICLE{WendtSPM2007,
%       author = {Herwig Wendt and Patrice Abry and St??phane Jaffard},
%       title = {Bootstrap for Empirical Multifractal Analysis},
%       journal = {IEEE Signal Proc. Mag.},
%       volume = {24},
%       number = {4},
%       pages = {38--48},
%       year = {2007},
%     }
%
% and (or)
%
%     @ARTICLE{WendtSP2009,
%       author = {Herwig Wendt and St??phane G. Roux and Patrice Abry and St??phane Jaffard},
%       title = {Wavelet leaders and bootstrap for multifractal analysis of images},
%       journal = {Signal Proces.},
%       volume = {89},
%       pages = {1100--1114},
%       year = {2009},
%     }
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(paramBS.Method)|(paramBS.Nresamp1<=1); doBS=0; else doBS=1; end
%% Parameter Setup
if isstruct(data) % data is structure functions
    doLS=0;
    logstatres=data;
    %--- which analysis are done
    try logstatres.DWT; if length(find(methodWT==1)); globalDWT=1;  else globalDWT=0;  end; catch globalDWT=0; methodWT=methodWT(find(methodWT~=1)); end
    try logstatres.LWT; if length(find(methodWT==2)); globalLWT=1;  else globalLWT=0;  end; catch globalLWT=0; methodWT=methodWT(find(methodWT~=2)); end
    % read parameters of LogScale
    j1min=logstatres.j1min;
    j2max=logstatres.j2max;
    Cum=logstatres.Cum;
    gamint=logstatres.gamint;
    MomNul=logstatres.MomNul;
    sym=logstatres.sym;
    EstFun=logstatres.EstFun;
    q=logstatres.q;
    Jflag=logstatres.Jflag;
    paramBS.Jflag=Jflag;
    T_S=logstatres.T_S;
    paramBS.T_S=T_S;
    j1=paramEST.j1; if length(j1)>1; j1=j1(1); end
    j2=paramEST.j2; if length(j2)>1; j2=j2(1); end
    if globalDWT;  logstatres.DWT.j1=round(j1);  logstatres.DWT.j2=min(round(j2),floor(j2max)); end;
    if globalLWT;  logstatres.LWT.j1=round(j1);  logstatres.LWT.j2=min(round(j2),floor(j2max)); end;
    scale1=exp(j1*log(2)); scale2=exp(min(j2,j2max)*log(2));
else
    doLS=1;
    %--- which analysis are to be done
    if length(find(methodWT==1)); globalDWT=1;  else globalDWT=0;  end
    if length(find(methodWT==2)); globalLWT=1;  else globalLWT=0; end
    %--- Estimation parameters
    j1=paramEST.j1; if length(j1)>1; j1=j1(1); end
    j2=paramEST.j2; if length(j2)>1; j2=j2(1); end
    MomNul=paramEST.MomNul;
    sym=paramEST.sym;
    Cum=paramEST.Cum;
    gamint=paramEST.gamint;
    EstFun=paramEST.Fun;
    q=paramEST.q;
    %--- BS Methods/Parameters
    Jflag=paramBS.Jflag;
    T_S=paramBS.T_S;
    if ~Jflag|~doBS; j1min=1; else; j1min=j1; end
end
% Nm=length(methodWT);  % number of analysing methods
if globalDWT&globalLWT; 
    Nm=2; methodWT=[1 2];
else
    if globalDWT; Nm=1; methodWT=[1]; end
    if globalLWT; Nm=1; methodWT=[2]; end
end
CI=paramBS.CI;
TEST=paramBS.TEST;
Alpha=paramBS.Alpha(1);
Type=paramTest.Type;
Tnull=paramTest.Tnull;
wtype=paramEST.wtype;
[paramBS] = correctparamBS(paramBS);

%--- Check what is to estimate
Fun=bin2dec(num2str(EstFun));
globalEstFun=EstFun;
globalFun=Fun;
if Fun>=4; ZQ=1; globalZQ=1; else; ZQ=0; globalZQ=0; end
if (Fun~=1) && (Fun~=4) && (Fun~=5); DH=1; globalDH=1; else;  DH=0; globalDH=0;end
if rem(Fun,2)==1; CP=1; globalCP=1; else; CP=0; globalCP=0; end

%--- Fig Numbers
%% --- ! --- param.EstFun
if FigNum>0
    FN=FigNum+[2:7];
    if ((verbose==2)|(verbose==3)|(verbose==21))&(ZQ|CP); PlotLogScale=1; else PlotLogScale=0; end
    if ((verbose==1)|(verbose==2)|(verbose==3)|(verbose==11)|(verbose==21)); PlotResult=1; else PlotResult=0; end
    figs=PlotLogScale+PlotResult;
else
    FN=zeros(1,5);
    figs=0;
    PlotLogScale=0; 
    PlotResult=0;
end
% if (FN>0)&((verbose==2)|(verbose==3)|(verbose==21))&(ZQ|CP); PlotLogScale=1; else PlotLogScale=0; end
% if (FN>0)&((verbose==1)|(verbose==2)|(verbose==3)|(verbose==11)|(verbose==21)); PlotResult=1; else PlotResult=0; end
if verbose>0
    if (verbose==1)|(verbose==2)|(verbose==3); TextResult=[1 CI TEST]; else; TextResult=[0 0 0]; end
    if verbose==3; doInter=1; else; doInter=0; end
else
    TextResult=[0 0 0];
    doInter=0;
end
%% Interactive mode parameters
INTERACT=1;
globalMethod=methodWT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- GET STRUCTURE FUNCTIONS
if doLS
    % check size of data
    [a,b]=size(data);
    if (a>1)&(b>1); imagedata=1; else; imagedata=0; end
    if ~imagedata
        if a>b; data=data'; end
%     else
%         if (a~=b)|(round(log2(a))~=log2(a)); error('Image must be square and its size a multiple of 2 (RICE toolbox)'); end
    end
    if globalDWT|globalLWT
        if verbose; disp('Discrete WT Analysis ...'); end;
        paramEST.j1=round(paramEST.j1); paramEST.j2=round(paramEST.j2);
        if ~imagedata
            [logstatDyad, logstatres.parBSestDyad] = MFA_BS_statlog1d(data, methodWT, paramEST, paramBS, Jflag, T_S);
        else
            [logstatDyad, logstatres.parBSestDyad] = MFA_BS_statlog2d(data, methodWT, paramEST, paramBS, Jflag, T_S);
        end
    end
    try logstatDyad.DWT; logstatres.DWT=logstatDyad.DWT; catch; end
    try logstatDyad.LWT; logstatres.LWT=logstatDyad.LWT; catch; end
    % get maximum possible scale
    if globalDWT; logstatres.DWT.imagedata=imagedata; end
    if globalLWT; logstatres.LWT.imagedata=imagedata; end
    j2max=[];
    try logstatres.DWT.scale; j2max=[j2max log2(logstatres.DWT.scale(end))]; catch; end;
    try logstatres.LWT.scale; j2max=[j2max log2(logstatres.LWT.scale(end))]; catch; end;
    j2max=min(j2max); j2=min(j2,j2max);
    % write parameters for LogScale
    logstatres.j1min=j1min;
    logstatres.j2max=j2max;
    logstatres.Cum=Cum;
    logstatres.gamint=gamint;
    logstatres.MomNul=MomNul;
    logstatres.sym=sym;
    logstatres.EstFun=EstFun;
    logstatres.q=q;
    logstatres.Jflag=Jflag;
    logstatres.T_S=T_S;
end
if ~isempty(find(methodWT==1))|~isempty(find(methodWT==2)); parBSestDyad=logstatres.parBSestDyad; end

while INTERACT  % interactive loop

    %--- GET FINAL ESTIMATES: LINEAR REGRESSIONS
    parBSestDyad.Alpha=Alpha;

    % -- j1, j2 may be non integer
    if globalDWT        
        [est.DWT]=MFA_BS_regrest(logstatres.DWT, parBSestDyad, paramEST, wtype);
    end
    if globalLWT
        [est.LWT]=MFA_BS_regrest(logstatres.LWT, parBSestDyad, paramEST, wtype);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CALCULATE CONFIDENCE LIMITS and TESTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--- DWT: CONFIDENCE LIMITS and TESTS ------------------------
    if doBS % ONLY IF BOOTSTRAP
        if globalDWT
            Nest=length(est.DWT.LEst);
            LE=[0 cumsum(est.DWT.LEst)];       % Length of estimates
            if globalFun>1; indZ=1:LE(end-1); indC=LE(end-1)+1:LE(end); else; indC=1:Cum; end %LE(5)=Cum;
            if CI
            for estid=1:LE(end)   % Loop through estimates: Cumulants 1 to Cum
                %%-- Prepare estimates for BS_CI.m (only 1 est at a time)
                [estW] = singleBSestimate(est.DWT, estid, parBSestDyad);
                %%-- Calculate Confidence Limits
                [confidence.DWT{estid}]=BS_CI(estW, parBSestDyad, Alpha);
            end
            end
            if TEST
            for estid=1:Cum   % Loop through estimates: Cumulants 1 to Cum
                %%-- Prepare estimates for BS_HT.m (only 1 est at a time)
                [estW] = singleBSestimate(est.DWT, indC(estid), parBSestDyad);
                %%-- Calculate Tests for Cumulants
                paramTest= struct('type', Type, 'T0', Tnull(estid));
                    [significance.DWT{estid}]=BS_HT(estW, parBSestDyad, paramTest, Alpha);
            end % Loop through estimates
            end
        end
        %--- LWT: CONFIDENCE LIMITS and TESTS ------------------------
        if globalLWT
            Nest=length(est.LWT.LEst);
            LE=[0 cumsum(est.LWT.LEst)];       % Length of estimates
            if globalFun>1; indZ=1:LE(end-1); indC=LE(end-1)+1:LE(end); else; indC=1:Cum; end
            if CI
            for estid=1:LE(end)   % Loop through estimates: Cumulants 1 to Cum
                %%-- Prepare estimates for BS_CI.m (only 1 est at a time)
                [estL] = singleBSestimate(est.LWT, estid, parBSestDyad);
                %%-- Calculate Confidence Limits
                    [confidence.LWT{estid}]=BS_CI(estL, parBSestDyad, Alpha);
            end
            end
            if TEST
            for estid=1:Cum   % Loop through estimates: Cumulants 1 to Cum
                %%-- Prepare estimates for BS_HT.m (only 1 est at a time)
                [estL] = singleBSestimate(est.LWT, indC(estid), parBSestDyad);
                %%-- Calculate Tests for Cumulants
                paramTest= struct('type', Type, 'T0', Tnull(estid));
                    [significance.LWT{estid}]=BS_HT(estL, parBSestDyad, paramTest, Alpha);
            end % Loop through estimates
            end
        end
    else
        if globalDWT; TEMP=est.DWT.LEst; elseif globalLWT; TEMP=est.LWT.LEst;    end;
        Nest=length(TEMP); LE=[0 cumsum(TEMP)];       % Length of estimates
        if globalFun>1; indZ=1:LE(end-1); indC=LE(end-1)+1:LE(end); else; indC=1:Cum;end
        confidence=NaN;
        significance=NaN;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% OUTPUT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if globalDWT|globalLWT
        varargout{1}=est;
        if CI; varargout{2}=confidence; else; varargout{2}=NaN; end
        if TEST; varargout{3}=significance; else; varargout{3}=NaN; end
        varargout{4}=logstatres;
    else
        varargout{1}=NaN;
        varargout{2}=NaN;
        varargout{3}=NaN;
        varargout{4}=NaN;
    end
    %return;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DISPLAY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constants
    bid1=max(1, floor(paramBS.Nresamp1*Alpha/2)); bid2=paramBS.Nresamp1+1-bid1;
    loge=log2(exp(1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FIGURES LOGSCALE
    %----- FIGURE LOGSCALE

    %% --- ! --- Input param.EstFun

    FontSize=12;
    FontSize1=8;
    if PlotLogScale
        % H min
        figure(FN(end)+1); clf; 
        if     methodWT(1)==1; logstattemp=logstatres.DWT; esttemp=est.DWT;  Mstr=['DWT'];
            elseif methodWT(1)==2; logstattemp=logstatres.LWT; esttemp=est.LWT;  Mstr=['LWT'];
            end;
            scales=log2(logstattemp.scale);
            % check scales that own bootstrap resamples
            scaleid=(find(scales>=round(j1min)));
            % check out regression range
            J1id=(logstattemp.j1); J2id=(logstattemp.j2);
            % for display
            jj1=logstattemp.j1; jj2=logstattemp.j2;
            
            sb=subplot(141); hold on;
            temp=log2(logstattemp.supcoefnoint);
            tempR=esttemp.h_minnoint; tempaest=esttemp.h_minnoint_aest;
            if tempR<=0; posX=0.7 ; else; posX=0.03; end
            text('string',[num2str(tempR)], 'fontsize',16,'FontName', 'Helvetica', 'Color', 'red','units','normalized', 'pos',[posX 0.9]);
                    text('string',['wtype ',num2str(wtype)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[posX 0.76]);
                    text('string',['j:[',num2str(jj1),'-',num2str(jj2),']'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[posX 0.69]);
                    grid on; plot(scales, temp, 'k.-'); plot(scales, scales*tempR+tempaest, 'b--'); 
                    plot(scales(J1id:J2id), scales(J1id:J2id)*tempR+tempaest, 'b'); hold off;
                    title([' h_{min} Coeff LogScale before fractional integration (\gamma=',num2str(logstattemp.gamint),')']); xlabel('j');
                    ylabel('log_2 sup [d(j,.)]');
             sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)');
            hold off;
            
             sb=subplot(142); hold on;
            temp=log2(logstattemp.supcoef);
            tempR=esttemp.h_min; tempaest=esttemp.h_min_aest;
            if tempR<=0; posX=0.7 ; else; posX=0.03; end
            text('string',[num2str(tempR)], 'fontsize',16,'FontName', 'Helvetica', 'Color', 'red','units','normalized', 'pos',[posX 0.9]);
                    text('string',['wtype ',num2str(wtype)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[posX 0.76]);
                    text('string',['j:[',num2str(jj1),'-',num2str(jj2),']'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[posX 0.69]);
                    grid on; plot(scales, temp, 'k.-'); plot(scales, scales*tempR+tempaest, 'b--'); 
                    plot(scales(J1id:J2id), scales(J1id:J2id)*tempR+tempaest, 'b'); hold off;
                    title([' h_{min} Coeff LogScale']); xlabel('j');
                    ylabel('log_2 sup [d(j,.)]');
             sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)');
            hold off;
            
            sb=subplot(143); hold on;
            temp=log2(logstattemp.supcoefL);
            tempR=esttemp.h_minL; tempaest=esttemp.h_minL_aest;
            if tempR<=0; posX=0.7 ; else; posX=0.03; end
            text('string',[num2str(tempR)], 'fontsize',16,'FontName', 'Helvetica', 'Color', 'red','units','normalized', 'pos',[posX 0.9]);
                    text('string',['wtype ',num2str(wtype)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[posX 0.76]);
                    text('string',['j:[',num2str(jj1),'-',num2str(jj2),']'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[posX 0.69]);
                    grid on; plot(scales, temp, 'k.-'); plot(scales, scales*tempR+tempaest, 'b--'); 
                    plot(scales(J1id:J2id), scales(J1id:J2id)*tempR+tempaest, 'b'); hold off;
                    title([' h_{min} Leaders LogScale']); xlabel('j');
                    ylabel('log_2 sup [L(j,.)]');
             sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)');
            hold off;
            
             sb=subplot(144); hold on;
            temp=log2(logstattemp.mincoef);
            tempR=esttemp.h_max; tempaest=esttemp.h_max_aest;
            if tempR<=0; posX=0.7 ; else; posX=0.03; end
            text('string',[num2str(tempR)], 'fontsize',16,'FontName', 'Helvetica', 'Color', 'red','units','normalized', 'pos',[posX 0.9]);
                    text('string',['wtype ',num2str(wtype)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[posX 0.76]);
                    text('string',['j:[',num2str(jj1),'-',num2str(jj2),']'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[posX 0.69]);
                    grid on; plot(scales, temp, 'k.-'); plot(scales, scales*tempR+tempaest, 'b--'); 
                    plot(scales(J1id:J2id), scales(J1id:J2id)*tempR+tempaest, 'b'); hold off;
                    title([' h_{max} Leaders LogScale']); xlabel('j');
                    ylabel('log_2 inf [L(j,.)]');
             sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)');
                    hold off;
                    
        figure(FN(1)); clf;
        for Mid=1:Nm
            if     methodWT(Mid)==1; logstattemp=logstatres.DWT; esttemp=est.DWT;  Mstr=['DWT'];
            elseif methodWT(Mid)==2; logstattemp=logstatres.LWT; esttemp=est.LWT;  Mstr=['LWT'];
            end;
            % check available scales
            scales=log2(logstattemp.scale);
            % check scales that own bootstrap resamples
            scaleid=(find(scales>=round(j1min)));
            % check out regression range
            J1id=(logstattemp.j1); J2id=(logstattemp.j2);
            % for display
            jj1=logstattemp.j1; jj2=logstattemp.j2;
            %Ltmp=cumsum([0 esttemp.LEst]);
            % Cumulants
            if CP
            figure(FN(1));
            for ii=1:Cum
                sb=subplot(Cum,Nm,Mid+(ii-1)*(Nm));estid=indC(ii);
                temp=logstattemp.est(estid,:)*loge;

                if doBS; tempB=sort(logstattemp.estB(:,estid,:))*loge; end;

                tempR=esttemp.t(estid); tempaest=esttemp.aest(estid); try esttemp.Q(estid); tempQ=esttemp.Q(estid); catch; tempQ=NaN; end;
                text('string',[num2str(tempR)], 'fontsize',FontSize,'FontName', 'Helvetica', 'Color', 'red','units','normalized', 'pos',[0.03 0.9]);
                text('string',['Q=',num2str(tempQ)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.83]);hold on;
                text('string',['wtype ',num2str(wtype)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.76]);
                text('string',['j:[',num2str(jj1),'-',num2str(jj2),']'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.69]);

                if doBS; errorbar(scales(scaleid), temp(scaleid), squeeze(tempB(bid1,:,(scaleid)))'-temp(scaleid), squeeze(tempB(bid2,:,(scaleid)))'-temp(scaleid), 'r'); end

                grid on; plot(scales, temp, 'k.-'); plot(scales, scales*tempR+tempaest, 'b--'); plot(scales(J1id:J2id), scales(J1id:J2id)*tempR+tempaest, 'b'); hold off;
                title([Mstr,' c_',num2str(ii),' LogScale']); xlabel('j');
                sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)');
            end
            end
            % Moments
            if ZQ
                figure(FN(Mid+1)); clf;
                lenq=length(q);
                sbid=zeros(1,3*ceil(lenq/3)); sbid(1:ceil(lenq/3))=1:3:3*ceil(lenq/3); sbid(ceil(lenq/3)+1:2*ceil(lenq/3))=2:3:3*ceil(lenq/3); sbid(2*ceil(lenq/3)+1:end)=3:3:3*ceil(lenq/3);
                for estid=1:lenq
                    sb=subplot(ceil(lenq/3),3,sbid(estid));
                    temp=logstattemp.est(estid,:);

                    if doBS; tempB=sort(logstattemp.estB(:,estid,:)); end;

                    tempR=esttemp.t(estid); tempaest=esttemp.aest(estid); try esttemp.Q(estid); tempQ=esttemp.Q(estid); catch; tempQ=NaN; end;
                    if q(estid)<=0; posX=0.7 ; else; posX=0.03; end
                    text('string',[num2str(tempR)], 'fontsize',FontSize,'FontName', 'Helvetica', 'Color', 'red','units','normalized', 'pos',[posX 0.9]);
                    text('string',['Q=',num2str(tempQ)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[posX 0.83]);hold on;
                    text('string',['wtype ',num2str(wtype)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[posX 0.76]);
                    text('string',['j:[',num2str(jj1),'-',num2str(jj2),']'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[posX 0.69]);

                    if doBS; errorbar(scales(scaleid), temp(scaleid), squeeze(tempB(bid1,:,(scaleid)))'-temp(scaleid), squeeze(tempB(bid2,:,(scaleid)))'-temp(scaleid), 'r'); end;

                    grid on; plot(scales, temp, 'k.-'); plot(scales, scales*tempR+tempaest, 'b--'); plot(scales(J1id:J2id), scales(J1id:J2id)*tempR+tempaest, 'b'); hold off;
                    title([Mstr,' \zeta(q) LogScale, q=',num2str(q(estid))]); xlabel('j');
                    sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)');
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FIGURES RESULTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% --- ! --- Input param.EstFun

    if PlotResult
        BoxW=0.5;
        LineWidth=0.5;
        MarkerSize=10;
        % Cumulants
        if CP
        hh1=figure(FigNum); clf;     figure(hh1);
        for Mid=1:Nm
            if     methodWT(Mid)==1; logstattemp=logstatres.DWT; esttemp=est.DWT;  Mstr=['DWT'];
            elseif methodWT(Mid)==2; logstattemp=logstatres.LWT; esttemp=est.LWT;  Mstr=['LWT'];
            end;
            % for display
            jj1=logstattemp.j1; jj2=logstattemp.j2;            
            for ii=1:Cum
                sb=subplot(Cum,Nm,Mid+(ii-1)*(Nm)); estid=indC(ii); cla;

                if doBS; boxplot([esttemp.T(:,(estid))], 'notch','on', 'widths', BoxW); hold on; end;

                plot([0.5 1-BoxW/2-0.05],[esttemp.t((estid)) esttemp.t((estid))], 'k', [1+BoxW/2+0.05 1.5],[esttemp.t((estid)) esttemp.t((estid))], 'k','LineWidth',LineWidth ); hold off;
                if doBS; miny=min(esttemp.T(:,(estid))); maxy=max(esttemp.T(:,(estid)));  else; miny=esttemp.t((estid))*0.9; maxy=esttemp.t((estid))*1.1; end
                dist=(maxy-miny)/8; miny=miny-dist; maxy=maxy+dist;
                title([Mstr,' c_',num2str(ii)]); %xlabel(['c_',num2str(estid)]); ylabel(''); grid on;
                xlabel([' ']); ylabel(''); grid on;
                set(gca, 'XTickLabel', ' ');
                text('string',[num2str(esttemp.t((estid)))], 'fontsize',FontSize,'FontName', 'Helvetica', 'Color', 'red','units','normalized', 'pos',[0.03 0.9]);
                text('string',['Frac Int ',num2str(gamint)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.83]);
                text('string',['Mnul ',num2str(MomNul)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.76]);
                text('string',['j:[',num2str(jj1),'-',num2str(jj2),']'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.69]);
                text('string',['wtype ',num2str(wtype)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.62]);
                sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)');
            end
            end
        end
        % Moments
        if ZQ|DH
            hh2=figure(FigNum+1); clf;
            figure(hh2);
            for Mid=1:Nm
                if     methodWT(Mid)==1; logstattemp=logstatres.DWT; esttemp=est.DWT;  Mstr=['DWT'];
                elseif methodWT(Mid)==2; logstattemp=logstatres.LWT; esttemp=est.LWT;  Mstr=['LWT'];
                end;
                % for display
                jj1=logstattemp.j1; jj2=logstattemp.j2;
                % Zeta(q)
                if ZQ
                sb=subplot(2,Nm,Mid);
                ZL=esttemp.t(LE(1)+1:LE(2));
                if doBS;
                    ZBL=sort(esttemp.T(:,LE(1)+1:LE(2)))';
                    errorbar(q, ZL, ZBL(:,bid1)'-ZL, ZBL(:,bid2)'-ZL, 'r-','LineWidth',LineWidth, 'MarkerSize', MarkerSize);  hold on;
                end;

                plot(q, ZL, 'k.-','LineWidth',LineWidth, 'MarkerSize', MarkerSize); grid on; hold on;
                xlabel('q'); ylabel('\zeta(q)'); title([Mstr,' \zeta(q)']);
                text('string',['Frac Int ',num2str(gamint)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.93]);
                text('string',['Mnul ',num2str(MomNul)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.86]);
                text('string',['j:[',num2str(jj1),'-',num2str(jj2),']'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.79]);
                text('string',['wtype ',num2str(wtype)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.72]);
                sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)');
                end
                if DH
                sb=subplot(2,Nm,Mid+Nm);
                grid on; hold on;
                DqL=esttemp.t(LE(1+ZQ)+1:LE(2+ZQ));
                hqL=esttemp.t(LE(2+ZQ)+1:LE(3+ZQ));
                if doBS;
                    DqBL=sort(esttemp.T(:,LE(1+ZQ)+1:LE(2+ZQ)))';
                    hqBL=sort(esttemp.T(:,LE(2+ZQ)+1:LE(3+ZQ)))';
                    DqLlo=DqBL(:,bid1)'; DqLhi=DqBL(:,bid2)';
                    hqLlo=hqBL(:,bid1)'; hqLhi=hqBL(:,bid2)';
                    for kq=2:length(q);
                        plot([hqL(kq) hqL(kq)],[DqLlo(kq) DqLhi(kq)],'r');
                        plot([hqLlo(kq) hqLhi(kq)],[DqL(kq) DqL(kq)],'r');
                    end;
                end
                plot(hqL(2:end), DqL(2:end), 'k.-','LineWidth',LineWidth, 'MarkerSize', MarkerSize); % Struct
                %axis([0 2 -0.5 1.5]); % flt
                hold off;
                xlabel('h'); title([Mstr,' D(h)']);
                text('string',['Frac Int ',num2str(gamint)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.93]);
                text('string',['Mnul ',num2str(MomNul)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.86]);
                text('string',['j:[',num2str(jj1),'-',num2str(jj2),']'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.79]);
                text('string',['wtype ',num2str(wtype)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.72]);
                sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)');
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DISPLAY TEXT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% --- ! --- Input param.EstFun

    if find(TextResult)
        Method=paramBS.Method;
        SigMeth=cell(1,8);
        SigMeth{1}='NOR ';
        SigMeth{2}='BAS ';
        SigMeth{3}='PER ';
        SigMeth{4}='STU ';
        SigMeth{5}='ADJB';
        SigMeth{6}='ADJP';
        TMeth=cell(1,4);
        TMeth{1}=' A: t>tnull';
        TMeth{2}=' A: t<tnull';
        TMeth{3}='|t-tnull|=0';
        TMeth{4}=' t-tnull=0 ';
        if TextResult(1)
            %% Estimates
            disp('********** * * * * * * * * * * * * * * * * * * * * * * * * * * * * **********');
            for Mid=1:Nm
                if     methodWT(Mid)==1; esttemp=est.DWT; Mstr=['DWT'];
                elseif methodWT(Mid)==2; esttemp=est.LWT; Mstr=['LWT'];
                end;
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                disp(['%%%%%     ESTIMATES ',Mstr,'      %%%%%']);
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                fprintf('[ j1 = %1g -- j2 = %2g -- wtype = %1g ]\n', j1, j2, wtype);
                disp('-----------------------------------');
                fprintf('  Fractional Integration: gamma=%1.2f \n', logstatres.gamint);
                fprintf('  Estimated h_min d_X before int : h_min=%1.2f \n', esttemp.h_minnoint);
                fprintf('  Estimated h_min d_X            : h_min=%1.2f \n', esttemp.h_min);
                fprintf('  Estimated h_min L_X            : h_min=%1.2f \n', esttemp.h_minL);
                fprintf('  Estimated h_max L_X            : h_max=%1.2f \n', esttemp.h_max);                
                disp('-----------------------------------');
                if doBS; fprintf('  Estimate                ( std* )\n'); else; fprintf('  Estimate  \n'); end
                if CP
                disp('-----------------------------------');
                for estid=1:Cum   % Loop through estimates: Cumulants 1 to Cum
                    if doBS
                        fprintf('        Cum %1g = %2.3f    (%2.3f)\n',estid,esttemp.t(indC(estid)),esttemp.stdt(indC(estid)));
                    else
                        fprintf('        Cum %1g = %2.3f  \n',estid,esttemp.t(indC(estid)));
                    end
                end                
                disp('--------------');
                end
                if ZQ
                    disp('-----------------------------------');
                    for estid=LE(1)+1:LE(2)   % Loop through estimates: Cumulants 1 to Cum
                        if doBS
                            fprintf('         zeta = %2.3f    (%2.3f)       q=%1g\n',esttemp.t(estid),esttemp.stdt(estid),q(estid));
                        else
                            fprintf('         zeta = %2.3f       q=%1g\n',esttemp.t(estid),q(estid));
                        end
                    end
                    disp('--------------');
                end
                if DH
                    disp('-----------------------------------');
                    for estid=LE(1+globalZQ)+1:LE(2+globalZQ)   % Loop through estimates: Cumulants 1 to Cum
                        if doBS;
                            fprintf('            D = %2.3f    (%2.3f)       q=%1g\n',esttemp.t(estid),esttemp.stdt(estid),q(estid-LE(1+globalZQ)));
                        else
                            fprintf('            D = %2.3f       q=%1g\n',esttemp.t(estid),q(estid-LE(1+globalZQ)));
                        end
                    end
                    disp('--------------');
                    disp('-----------------------------------');
                    for estid=LE(2+globalZQ)+1:LE(3+globalZQ)   % Loop through estimates: Cumulants 1 to Cum
                        if doBS
                            fprintf('            h = %2.3f    (%2.3f)       q=%1g\n',esttemp.t(estid),esttemp.stdt(estid),q(estid-LE(2+globalZQ)));
                        else
                            fprintf('            h = %2.3f       q=%1g\n',esttemp.t(estid),q(estid-LE(2+globalZQ)));
                        end
                    end
                    disp('--------------');
                end
                disp('  ');
            end
        end % if Text
        if TextResult(2)&doBS&CI
            %% Confidence
            disp('********** * * * * * * * * * * * * * * * * * * * * * * * * * * * * **********');
            for Mid=1:Nm
                if     methodWT(Mid)==1; esttemp=est.DWT; conftemp=confidence.DWT;  Mstr=['DWT'];
                elseif methodWT(Mid)==2; esttemp=est.LWT; conftemp=confidence.LWT;  Mstr=['LWT'];
                end;
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                disp(['%%%%%      CONFIDENCE INTERVALS ',Mstr,'     %%%%%']);
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                    for methodid=1:length(Method)
                        MID=Method(methodid);
                        disp('  ');
                        fprintf('---  METHOD=%s   - ALPHA=%1.2f     ---\n',SigMeth{MID}, Alpha);
                        fprintf('  Estimate       [CIlo    CIhi]  (estimate (std*))\n');
                        if CP
                            disp('-----------------------------------');
                            for estid=1:Cum   % Loop through estimates: Cumulants 1 to Cum
                                fprintf('   Cumulant %1g : [%2.3f   %2.3f]  ( %2.3f (%2.3f) )\n',estid, conftemp{indC(estid)}{MID}.lo, conftemp{indC(estid)}{MID}.hi, esttemp.t(indC(estid)),esttemp.stdt(indC(estid)));
                            end
                            disp('--------------');
                        end
                        if ZQ
                            disp('-----------------------------------');
                            for estid=LE(1)+1:LE(2)   % Loop through estimates: Cumulants 1 to Cum
                                fprintf('         zeta : [%2.3f   %2.3f]  ( %2.3f (%2.3f) )   q=%1g\n', conftemp{(estid)}{MID}.lo, conftemp{(estid)}{MID}.hi, esttemp.t((estid)),esttemp.stdt((estid)),q(estid));
                            end
                            disp('--------------');
                        end
                        if DH
                            disp('-----------------------------------');
                            for estid=LE(1+globalZQ)+1:LE(2+globalZQ)   % Loop through estimates: Cumulants 1 to Cum
                                fprintf('             D : [%2.3f   %2.3f]  ( %2.3f (%2.3f) )   q=%1g\n', conftemp{(estid)}{MID}.lo, conftemp{(estid)}{MID}.hi, esttemp.t((estid)),esttemp.stdt((estid)),q(estid-LE(1+globalZQ)));
                            end
                            disp('--------------');
                            disp('-----------------------------------');
                            for estid=LE(2+globalZQ)+1:LE(3+globalZQ)   % Loop through estimates: Cumulants 1 to Cum
                                fprintf('             h : [%2.3f   %2.3f]  ( %2.3f (%2.3f) )   q=%1g\n', conftemp{(estid)}{MID}.lo, conftemp{(estid)}{MID}.hi, esttemp.t((estid)),esttemp.stdt((estid)),q(estid-LE(2+globalZQ)));
                            end
                            disp('--------------');
                        end
                    end
                    disp('  ');
            end
        end % if Text
        if TextResult(3)&doBS&TEST
            %% Significance
            if CP
            disp('********** * * * * * * * * * * * * * * * * * * * * * * * * * * * * **********');
            for Mid=1:Nm
                if     methodWT(Mid)==1; esttemp=est.DWT; conftemp=confidence.DWT; signiftemp=significance.DWT; Mstr=['DWT'];
                elseif methodWT(Mid)==2; esttemp=est.LWT; conftemp=confidence.LWT; signiftemp=significance.LWT;  Mstr=['LWT'];
                end;
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                disp(['%%%      CUMULANTS:        TESTS ',Mstr,'     %%%%']);
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                    for estid=1:Cum
                        disp(' ');
                        fprintf('----- Cumulant %1g   - ALPHA=%1.2f     ---\n', estid, Alpha);
                        fprintf('Hnull: c_%1g,0 =%2.4f  -- c_%1g = %2.4f  ( %2.4f )\n',estid,Tnull(estid),estid, esttemp.t(indC(estid)),esttemp.stdt(indC(estid)));
                        for typeid=1:length(Type)
                            TID=Type(typeid);
                            fprintf('- Test Statistic %s        \n',  TMeth{TID});
                            fprintf('Method     Reject   P-value     \n');
                            for methodid=1:length(Method)
                                MID=Method(methodid);
                                fprintf('%s :     %1g         %2.3f \n',SigMeth{MID}, signiftemp{estid}{MID}{TID}.reject, signiftemp{estid}{MID}{TID}.plevel );
                            end
                        end
                    end
            end
            end
        end % if Text
        disp('********** * * * * * * * * * * * * * * * * * * * * * * * * * * * * **********');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Interaction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% --- ! --- Input param.EstFun

    if doInter
        %--- cp only or as well zeta(q) and D(h)
        CHANGE=0;
        disp('  ');
        disp('----------------------------------------------');
        disp('            INTERACTIVE MODE ');
        disp('----------------------------------------------');
        disp('  ');
        if globalZQ
            temp=input('** Estimate zeta(q) ? [0 1] :   ');
            if ~isempty(temp);        
                oldZQ=ZQ;
                if temp(1)>0; ZQ=1; else; ZQ=0;  end; 
                if oldZQ~=ZQ; CHANGE=1; end
            end
        end
        if globalDH
            temp=input('** Estimate D(h) ? [0 1] :   ');
            if ~isempty(temp);        
                oldDH=DH;
                if temp(1)>0; DH=1; else; DH=0;  end; 
                if oldDH~=DH; CHANGE=1; end
            end            
        end
        if globalCP
            temp=input('** Estimate cp ? [0 1] :   ');
            if ~isempty(temp);        
                oldCP=CP;
                if temp(1)>0; CP=1; else; CP=0;  end; 
                if oldCP~=CP; CHANGE=1; end
            end            
        end
        % j1, j2, wtype
        temp=input(['** New [j1 j2 (wtype)] (j1>=',num2str(j1min),' j2<=',num2str(j2max),') ? ( Current [',num2str(max(j1,j1min)),' ',num2str(min(j2,j2max)),' ',num2str(wtype),']) :   ']);
        if  ~isempty(temp) ;
            CHANGE=1;
            try temp(1); j1=max(j1min, temp(1)); catch; end;
            try temp(2); j2=min(j2max, temp(2)); catch; end;
            try temp(3); if(temp(3)==0)|(temp(3)==1)|(temp(3)==2); wtype=temp(3);end;  catch; end;
            if globalDWT;  logstatres.DWT.j1=round(j1);  logstatres.DWT.j2=round(j2); end;
            if globalLWT;  logstatres.LWT.j1=round(j1);  logstatres.LWT.j2=round(j2); end;
        end
        % Analysis Method
        temp=input(['** Analysis Methods [',num2str(globalMethod),'] ? :   ']);
        if ~isempty(temp)
            CHANGE=1;
            methodWT=[];
            for i=1:length(globalMethod)
                if find(temp==globalMethod(i)); methodWT=[methodWT globalMethod(i)]; end
            end
            if isempty(methodWT); methodWT=globalMethod; disp(['Methods not valid. Keeping [',num2str(globalMethod),'] ...']); end
        end
        Nm=length(methodWT);
        % Alpha
        if doBS;
            temp=input(['** New alpha ? (current ',num2str(Alpha),') :   ']);
            if ~isempty(temp)
                CHANGE=1;
                if (temp(1)>0)&(temp(1)<1); Alpha=temp(1); else; disp(['New alpha not valid, keeping ',num2str(Alpha),' ...']); end
            end
        end
        % Text Display
        temp=input(['** Display tables [Estimates Confidence Tests] ? Current [',num2str(TextResult),'] :   ']);
        if ~isempty(temp)
            CHANGE=1;
            for i=1:3
                try temp(i); if (temp(i)==0)|(temp(i)==1); TextResult(i)=temp(i); end; catch; end
            end
            if TextResult(2)==1; CI=1; end;
            if TextResult(3)==1; TEST=1; end;
        end
        % Figures
        disp(['** Display Figures (+ FigNum -optional)? [Figure (FigNum)] ? Current [',num2str(figs),' (',num2str(FigNum),')]']);
        if CHANGE==0; disp(' ---> Press ENTER to exit Interactive Mode <---'); end
        temp=input(['       Figure=( 0-off 1-estimates 2-LogScale)   :   ']);
        if isempty(temp)
            if ~CHANGE
                disp('EXITING .');
            end
        else
            CHANGE=1;
            try temp(1); if (temp(1)==0)|(temp(1)==1)|(temp(1)==2); figs=temp(1); end; catch; end
            try temp(2); if (temp(2)>=1); FigNum=round(temp(2)); end; catch; end
            if (figs>0)&(FigNum==0); FigNum=1; end
        end
        if FigNum>0
            FN=FigNum+[2:7];
        end
        if (FigNum>0)&(verbose>1); PlotLogScale=1; else PlotLogScale=0; end
        if (FigNum>0)&(verbose>0); PlotResult=1; else PlotResult=0; end
        if ~CHANGE; INTERACT=0; end
    end
    if ~doInter; INTERACT=0; end
end







