%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     DEMO OF BOOTSTRAP MFA ANALYSIS with MF_BS_test_constancy.m:
%                demo_MF_BS_test_constancy.m
%     Wavelet domain double block bootstrap based test of hypothesis of
%     time constancy of multifractal attributes
%
% Herwig Wendt, Lyon, 2006 - 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Please refer to the references on my web site (and the below given
% references) for details on the method.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this tool, please quote:
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
% and:
%
%     @INPROCEEDINGS{WendtICASSP2008,
%       author = {Herwig Wendt and Patrice Abry},
%       title = {Bootstrap tests for the time constancy of multifractal attributes},
%       booktitle = {Proc. IEEE Int. Conf. Acoust., Speech, and Signal Proc. (ICASSP)},
%       address = {Las Vegas, USA},
%       year = {2008},
%     }
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
pathMF_BS_toolbox;         % add subpaths of toolbox

%% CHOOSE EXAMPLE
H0=1; 
% 1: data (mrw) where multifractal attributes are constant with time
% 0: data (mrw) where multifractal attributes are not constant with time:
%       concatenation of two realizations of mrw with different parameters


if H0; 
    load H0_mrw_075_008_n8192; H0=1; TSTR=['mrw with c_p=[0.75,-0.08] - constant'];
else
    load H1_mrw_075068_008001_n8192; H0=0;  TSTR=['concatenation of mrw with c_p=[0.75,-0.08] and c_p=[0.68,-0.01] - not constant'];
end
figure(11111); clf; plot([data]); grid on; title(TSTR);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  PARAMETER SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- Number of time windows to be used
M=2;                % (note: trade off power <--> time resolution)

%--- Analysis Methods to use
methodWT=[2];       % 1 - DWT  (discrete wavelet transform coefficients)
                    % 2 - LWT  (wavelet Leaders)
%--- Estimation Parameters
MomNul=3;           % vanishing moments wavelet (Daubechies')
Norm=1;             % Normalization of wavelet coefficients 
j1=2;               % j1,j2 - scaling range
j2=log2(N)+1-log2(M)-5;
wtype=1;            % linear regression : 0 - ordinary ; 1 - weighted nj(j) ; 2 - weighted with estimated variance
Fun=001;            %       number xyz - what is to be calculated:  
                    %       x : zeta(q) [0 or 1]    (scaling exponents)
                    %       y : D(h) [0 or 1]       (multifractal spectrum)
                    %       z : cp [0 or 1]         (log-cumulants)
                    %       e.g. param.EstFun=001 calculates only Cp, 
                    %            param.EstFun=110 calculates SF zeta(q) and D(h)

Cum=3;              % Highest order of log cumulants cp to estimate
q=[1 2];           % Moments to estimate (zeta(q) and D(h))

%--- BOOTSTRAP
T_S=0;              % time - scale block [1] or time block [0] bootstrap ([1]  recommended)
% B1=29;              % # primary resamples
% B2=49;              % # double bootstrap resamples (double bootstrap is necessary for calculating the test statistic!)
B1=9;              % # primary resamples
B2=9;              % # double bootstrap resamples (double bootstrap is necessary for calculating the test statistic!)
if ~T_S             % Block Length
    Block=2*MomNul;     % Block Length
else
    Block=floor(N/32);
end
Method=[3];         % Methods to be used [1-Normal; 2-Basic; 3-Percentile; 4-Studentised; 5-Adjusted Basic; 6-Adjusted Percentile]
                    %      3 is recommended
Alpha=0.1;          % Significance Level (1-Alpha)
Jflag=1;            % 1: do not do Bootstrap on scales that are not implied in regression (use 0 only if you need confidence intervals for scales not involved in regressions)

%--- VERBOSITY PARAMETERS
FigNum=1;           % 0 - no Figures
verbose=1;          % Verbosity level:
                    %   0 : batch mode
                    %   1 : table of test decisions and p-values 
                    %   2 : table of test decisions , p-values, test statistic and critical value
                    %       Figure estimate if FigNum
                    %   3 : 2 + table with estimates (global and windowed)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rand('state',sum(100*clock));
% check and write parameters:
checkParam_MF_BS_test_constancy;
[paramBS, paramEst] = MF_BS_test_constancy_param(B1,B2,Block,Method,T_S,Alpha,Fun,j1,j2,wtype,Jflag,MomNul,methodWT,Cum,q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
[test, estm, logstatm] = MF_BS_test_constancy(data, M, paramEst, paramBS, verbose, FigNum);
toc;

% OUTPUT:
%
% estimates are organized as [zeta(q), D(q), h(q), cp]
%-- test: test.DWT, test.LWT
%%  test.DWT.t(estid)            - test statistics
%   test.DWT.T(bsid,estid)       - bootstrap test statistics
%   test.DWT.d(estid)            - test decisions
%   est.DWT.p(estid)             - test p-values
%
%%  estm.DWT.t(m,estid)          - estimates over windows
%   estm.DWT.t_all(estid)        - global estimates
%   estm.DWT.T(m,bsid,estid)     - bootstrap estimates (per window)
%   estm.DWT.stdt(m,estid)       - bootstrap standard deviation of estimtes (per window)
%   estm.DWT.stdT(m,bsid,estid)  - bootstrap estimates  (per window)
%   estm.DWT.nj(m,j)             - number of coefficients per scale (per window)
%   estm.DWT.njall(j)            - number of coefficients per scale (global)
%   estm.DWT.j1                  - scaling range [j1,j2] used (windows)
%   estm.DWT.j2
%   estm.DWT.j1all               - scaling range [j1,j2] used (global)
%   estm.DWT.j2all
%   estm.DWT.LEst                - length of estimates zeta(q), D(q), h(q),
%
%-- logstat
%%  logstat.DWT.Sqj(m,estid,j)          - Structure Functions
%   logstat.DWT.SqjBcut(m,bsid,estid,j) - Boostrap-then-cut Structure Functions (windows)
%   logstat.DWT.SqjBall(m,bsid,estid,j) - Cut-then-Boostrap Structure Functions (windows)
