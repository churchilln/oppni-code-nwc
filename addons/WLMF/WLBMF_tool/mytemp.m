%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     DEMO OF BOOTSTRAP MFA ANALYSIS with MF_BS_tool.m:
%                demo_MF_BS_tool.m
%     Wavelet coefficient and leader based estimation of multifractal attributes
%     Wavelet domain block bootstrap based confidence limits and significance tests
% 
% Herwig Wendt, Lyon, 2006 - 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Please refer to the references on my web site (and the below given
% references) for details on the method.
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

clear all;
pathMF_BS_toolbox;         % add subpaths of toolbox

%% CHOOSE EXAMPLE
MFprocess = 3 ;
% case 1: fbm,       N=4096         (1d, self-similar/monofractal)
% case 2: mrw,       N=32768        (1d, multifractal)
% case 3: cpc ln, N=1024 x 1024     (2d, multifractal)

addpath /Applications/MATLAB/work/WLBMF_example_data/
addpath /Applications/MATLAB/work/IRIT/ricedsp-rwt/src/

switch MFprocess
    case 1; % fbm, N=4096
        load fbm08n4096; N=length(data); figure(11111); clf; plot(data);  TSTR=['A realization of fractional Brownian motion (H=0.8, N=4096)'];
    case 2; % mrw, N=32768
        load mrw07005n32768;  N=length(data); figure(11111); clf; plot(data);  TSTR=['A realization of multifractal random walk (c_1=0.75, c_2=-0.05, N=32768)'];
    case 3; % cpc ln, N=1024x1024
        load cmcLN2d_00125_0025_n1024; N=length(data); figure(11111); clf; surf(data);  TSTR=['A realization of 2d compound poisson cascade with log-normal multipliers (c_1=0.0125, c_2=-0.025, N=1024x1024)'];
        data=data(1:768,3:999);
end
title(TSTR); grid on;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  PARAMETER SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- Analysis Methods to use
methodWT=[1 2];     % 1 - DWT  (discrete wavelet transform coefficients)
                    % 2 - LWT  (wavelet Leaders)
%--- Estimation Parameters
MomNul=3;           % vanishing moments wavelet (Daubechies')
gamint=0.7;           % Wavelet domain fractional integration by gamint
j1=3;               % j1,j2 - scaling range
j2=7;
wtype=1;            % linear regression : 0 - ordinary ; 1 - weighted nj(j) ; 2 - weighted with estimated variance
Fun=111;            %       number xyz - what is to be calculated:  
                    %       x : zeta(q) [0 or 1]    (scaling exponents)
                    %       y : D(h) [0 or 1]       (multifractal spectrum)
                    %       z : cp [0 or 1]         (log-cumulants)
                    %       e.g. param.EstFun=001 calculates only Cp, 
                    %            param.EstFun=110 calculates SF zeta(q) and D(h)

Cum=3;              % Highest order of log cumulants cp to estimate
q=[-10:-1 -0.5 0.5 1:10];           % Moments to estimate (zeta(q) and D(h))
    
%--- BOOTSTRAP
% Bootstrap parameters
T_S=0;              % time - scale block [1] or time block [0] bootstrap
B1=49;              % # primary resamples - NO BOOTSTRAP (ESTIMATION ONLY): B1=0;
B2=0;               % # double bootstrap resamples
CI=0;               % calculate confidence intervals
TEST=0;             % calculate tests
if T_S              % Block Length
    Block=floor(N/32);     
else
    Block=2*MomNul; 
end
Method=[3];         % Methods to be used [1-Normal; 2-Basic; 3-Percentile; 4-Studentised; 5-Adjusted Basic; 6-Adjusted Percentile]
Alpha=0.1;          % Significance Level (1-Alpha)
Jflag=1;            % 1: do not do Bootstrap on scales that are not implied in regression (use 0 only if you need confidence intervals for scales not involved in regressions)
% Test Parameters
Tnull=zeros(1,Cum); % Null Hypothesis Parameter (for cp only)
Type=[4];           % Null Hypothesis H0: T=Tnull Types [  HA: T>Tnull;   HA: T<Tnull;   H0: |T-Tnull|=0;   H0: T-Tnull=0  ]

%--- VERBOSITY PARAMETERS
FigNum=1;           % 0 - no Figures
verbose=2;          % Verbosity level:
                    %   0 : batch mode
                    %   1 : estimate, CI, Test table
                    %       Figure estimate if FigNum
                    %   11: Figure estimate if FigNum
                    %   2 : 1 + Figure LogScale if FigNum
                    %   21: 11+ Figure LogScale if FigNum
                    %   3 : 2 + interactive mode (change Estimates, j1, j2, wtype,
                    %           Analysis method, Alpha, Text output, Figures)

                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rand('state',sum(100*clock));
% check and write parameters:
checkParam_MF_BS_tool;
[paramEST, paramBS, paramTest]=MF_BS_tool_param(MomNul,gamint,j1,j2,wtype,Fun,Cum,q,B1,B2,Block,Method,Alpha,Jflag,Type,Tnull,T_S,CI,TEST);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
[est,conf,signif,logstat] = MF_BS_tool(data, paramEST, methodWT, paramBS, paramTest, verbose, FigNum);
toc;
% OUTPUT:
%
% estimates are organized as [zeta(q), D(q), h(q), cp]
%-- est: est.DWT, est.LWT
%%  est.DWT.t(estid)            - estimates
%   est.DWT.stdt(estid)         - bootstrap standard deviation
%   est.DWT.T(bsid,estid)       - bootstrap estimates
%   est.DWT.aest(estid)         - regression instersect estimation
%   est.DWT.Q(estid)            - regression quality estimation
%   est.DWT.LEst                - length of estimates zeta(q), D(q), h(q), cp
%
%-- conf: conf.DWT, conf.LWT
%   conf.DWT{estid}{Methodid}.name  - Method name
%%  conf.DWT{estid}{Methodid}.lo    - lower confidence limit
%%  conf.DWT{estid}{Methodid}.hi    - upper confidence limit
%
%-- signif: signif.DWT, signif.LWT
%   signif.DWT{estid}{Methodid}{Typeid}.name      - Method name
%   signif.DWT{estid}{Methodid}{Typeid}.reject    - Test decision
%%  signif.DWT{estid}{Methodid}{Typeid}.plevel    - plevel
%
%-- logstat
%%  logstat.DWT.est(estid, j)           - Structure Functions
%   logstat.DWT.nj(j)                   - number of coefficients at each scale
%   logstat.DWT.estB(bsid, estid, j)    - Bootstrap structure functions
%   logstat.DWT.Vest(estid, j)          - Bootstrap standard deviation
%   logstat.DWT.j1                      - scaling range [j1,j2] used
%   logstat.DWT.j2
%   logstat.DWT.scale                   - scale 2^j
%   logstat.DWT.imagedata               - imaghe or time series
%   logstat.parBSestDyad            - bootstrap parameters used
%   logstat.j1min                   - minimal scale for regression (bootstrap)
%   logstat.j2max                   - maximal scale for regression
%   logstat.Cum                     - number of cumulants calculated
%   logstat.Norm                    - Wavelet norm used
%   logstat.MomNul                  - vanishing moments wavelet used
%   logstat.EstFun                  - Estimation fuction parameter
%   logstat.q                       - moments q
%   logstat.Jflag                   - Jflag: Bootstrap on scale not used in regression
%   logstat.T_S                     - time scale block bootstrap
%
%% Minimum and maximum regularity estimation
% -- est.DWT.
% h_minnoint / h_min_aestnoint      - minimum Hoelder exponent before intergation (original image) and its regression intersect) 
% h_min / h_min_aest                - minimum Hoelder exponent after intergation
% h_minL / h_minL_aest              - minimum Hoelder exponent after intergation: Leaders (EXPERIMENTAL) 
% h_max / h_max_aest                - largest Hoelder exponent after intergation (Leaders based estimation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS bis 
% You can reuse structure functions stored in "LogScale" to test different
% hypothesis interactively
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('youp');
pause;

FigNum=101;         
verbose=1;          
Type=[4];     
%Example: Test a different Null Hypothesis
Tnull=Cp(1:Cum); % Null Hypothesis Parameter (for cp only)

% write parameters
checkParam_MF_BS_tool;
[paramEST, paramBS, paramTest]=MF_BS_tool_param(MomNul,Norm,j1,j2,wtype,Fun,Cum,q,B1,B2,Block,Method,Alpha,Jflag,Type,Tnull,T_S);
% ! only difference: "logstat" replaces "data"
tic;
[est2,conf2,signif2,logstat2] = MF_BS_tool(logstat, paramEST, methodWT, paramBS, paramTest, verbose, FigNum);
toc;
