function [out] = hurst_wavelet_mf( X, j1,j2, Qmax )
%
% Wavelet multifractal estimation -- uses method of Wendt et al.
%
% syntax:
%           Hurst = hurst_wavelet_mf(X)
%
% inputs:
%           X = data matrix, arranged (samples x time)
%
% outputs:
%           Cumul = matrix of dimensions (samples x 3), giving log-cumulants
%                   of orders 1-3 as column vectors [c1 c2 c3]. 
%                   c1 (column #1) is approximately equivalent to monofractal hurst exponent

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PARAMETER SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--- Analysis Methods to use
methodWT=[1 2];     % 1 - DWT  (discrete wavelet transform coefficients)
                    % 2 - LWT  (wavelet Leaders)

%--- Estimation Parameters
MomNul=3;           % vanishing moments wavelet (Daubechies')
gamint=0.99;        % Wavelet domain fractional integration by gamint (0.7 --> 0.9)
% j1=2;               % j1,j2 - scaling range [3,7]
% j2=4;

wtype=1;            % linear regression : 0 - ordinary ; 1 - weighted nj(j) ; 2 - weighted with estimated variance
Fun=111;            %       number xyz - what is to be calculated:  
                    %       x : zeta(q) [0 or 1]    (scaling exponents)
                    %       y : D(h) [0 or 1]       (multifractal spectrum)
                    %       z : cp [0 or 1]         (log-cumulants)
                    %       e.g. param.EstFun=001 calculates only Cp, 
                    %            param.EstFun=110 calculates SF zeta(q) and D(h)

Cum=3;              % Highest order of log cumulants cp to estimate
q=[-Qmax:-1 -0.5 0.5 1:Qmax];           % Moments to estimate (zeta(q) and D(h))
    
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
verbose=0;          % Verbosity level:
                    %   0 : batch mode
                    %   1 : estimate, CI, Test table
                    %       Figure estimate if FigNum
                    %   11: Figure estimate if FigNum
                    %   2 : 1 + Figure LogScale if FigNum
                    %   21: 11+ Figure LogScale if FigNum
                    %   3 : 2 + interactive mode (change Estimates, j1, j2, wtype,
                    %           Analysis method, Alpha, Text output, Figures)


checkParam_MF_BS_tool;

[paramEST, paramBS, paramTest]=MF_BS_tool_param(MomNul,gamint,j1,j2,wtype,Fun,Cum,q,B1,B2,Block,Method,Alpha,Jflag,Type,Tnull,T_S,CI,TEST);

for(i=1:size(X,1))
%     disp(['running ',num2str(i),'/',num2str(size(X,1))]);
    [hh1(i,:) DD1(i,:) Cumul(i,:)] = MF_BS_mod(X(i,:), paramEST, methodWT, paramBS, paramTest, verbose, FigNum);
end

out.Cumul = Cumul;
out.hh1   = hh1;
out.DD1   = DD1;

% figure;
% subplot(1,2,1), plot( mean(hh1,1)', mean(DD1,1)', '.-' );
% ylabel('singularity dimension D(h)');xlabel('singularity exponent h');
% title('multifractal spectrum');
% subplot(1,2,2), boxplot( Cumul ); 
% ylabel('cumulants');
% xlabel('moments');
