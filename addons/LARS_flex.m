function [beta, A, mu, C, c, gamma] = LARS_flex(X, Y, option, t, standardize, SparseCut)

% Least Angle Regression (LAR) algorithm.
% Ref: Efron et. al. (2004) Least angle regression. Annals of Statistics.
% option = 'lar' implements the vanilla LAR algorithm (default);
% option = 'lasso' solves the lasso path with a modified LAR algorithm.
% t -- a vector of increasing positive real numbers. If given, LARS stops and 
% returns the solution at t.
%
% Output:
% A -- a sequence of indices that indicate the order of variable inclusions;
% beta: history of estimated LARS coefficients;
% mu -- history of estimated mean vector;
% C -- history of maximal current absolute corrrelations;
% c -- history of current corrrelations;
% gamma: history of LARS step size.
% Note: history is traced by rows. If t is given, beta is just the
% estimated coefficient vector at the constraint ||beta||_1 = t.
%
% Remarks:
% 1. LARS is originally proposed to estimate a sparse coefficient vector in
% a noisy over-determined linear system. LARS outputs estimates for all
% shrinkage/constraint parameters (homotopy).
%
% 2. LARS is well suited for Basis Pursuit (BP) purpose in the real case. This lars function
% automatically terminates when the current correlations for inactive set are
% all zeros. The recovered coefficient vector is the last column of beta 
% with the *lasso* option. Hence, this function provides a fast and 
% efficient solution for the ell_1 minimization problem. 
% Ref: Donoho and Tsaig (2006). Fast solution of ell_1 norm minimization problems when the solution may be sparse.
% 
% Least angle regression (LAR) algorithm
% Author: Xiaohui Chen (xiaohuic@ece.ubc.ca)
% Version: 2012-Feb

if nargin < 5, standardize = true; disp('auto-stand'); end
if nargin < 4, t = Inf;            disp('no list!!!'); end
if nargin < 3, option = 'lar';     disp('LAR only!!'); end

if strcmpi(option, 'lasso'), lasso = 1; else, lasso = 0; end

eps = 1e-10;    % Effective zero

[n,p] = size(X); % n samples x p variables
if standardize,
    X = normalize(X);
    Y = Y-mean(Y);
end
m = min(p,n-1);     % Maximal number of variables in the final active set
T = length(t);      % number of stages recorded, defined by t-vector

beta = zeros(1,p);  % Beta regression coefficients
mu = zeros(n,1);    % Mean vector
gamma = [];         % LARS step lengths
A = [];             % declare "empty" list of active-set
Ac = 1:p;           % all variables currently inactive
nVars = 0;          % no variables in active set yet
signOK = 1;
i = 0;
mu_old = zeros(n,1);
t_prev = 0;
beta_t = zeros(T,p);
ii = 1;
tt = t;

% %% HACK: terminate if the cap on #desired variables reached
if( SparseCut > 0 ) maxVarNum = min( [SparseCut m] );
else                maxVarNum = m;
end

% LARS loop - add successive variables to the model
while( nVars < maxVarNum ) % HACK: originally nVars < m
    
    i = i+1;              % Step increment
    c = X'*(Y-mu);        % Current correlation - each var w/ residual (y-mu)
    C = max(abs(c));      % Maximal absolute correlation
    
    if C < eps || isempty(t), break; end    % Early stopping criteria
    if i == 1, addVar = find(abs(c)==C); end
    if signOK,            % * This condition only applies to LASSO model *
        A = [A,addVar];   % Add one variable to active set
        nVars = nVars+1;  % increment # variables
    end
    
    s_A = sign(c(A));     % sign of variables' correlation with y, in active set
    Ac = setdiff(1:p,A);  % Inactive set
    nZeros = length(Ac);  % declare number of non-active set
    X_A = X(:,A);         % matrix of active-set variables
    G_A = X_A'*X_A;       % Gram matrix (~covariance)
    invG_A = inv(G_A);    % inverse cov
    
    % scaling weight
    L_A = 1/sqrt(s_A'*invG_A*s_A);
    
    w_A = L_A*invG_A*s_A;   % Coefficients of equiangular vector u_A
    u_A = X_A*w_A;          % Equiangular vector
    a = X'*u_A;             % Angles between x_j and u_A
    
    beta_tmp = zeros(p,1);
    gammaTest = zeros(nZeros,2);
    if nVars == m,
        gamma(i) = C/L_A;   % Move to the least squares projection
    else
        for j = 1:nZeros,
            jj = Ac(j);
            gammaTest(j,:) = [(C-c(jj))/(L_A-a(jj)), (C+c(jj))/(L_A+a(jj))];
        end
        [gamma(i) min_i min_j] = minplus(gammaTest);
        addVar = unique(Ac(min_i));
    end
    beta_tmp(A) = beta(i,A)' + gamma(i)*w_A;    % Update coefficient estimates
    
    % Check the sign feasibility of lasso
    if lasso,
        signOK = 1;
        gammaTest = -beta(i,A)'./w_A;
        [gamma2 min_i min_j] = minplus(gammaTest);
        if gamma2 < gamma(i),   % The case when sign consistency gets violated
            gamma(i) = gamma2;
            beta_tmp(A) = beta(i,A)' + gamma(i)*w_A;    % Correct the coefficients
            beta_tmp(A(unique(min_i))) = 0;
            A(unique(min_i)) = [];  % Delete the zero-crossing variable (keep the ordering)
            nVars = nVars-1;
            signOK = 0;
        end
    end
    
    if  t(1) ~= Inf,
        t_now = norm(beta_tmp(A),1);

        if t_prev < t(1) && t_now >= t(1),
            beta_t(ii,A) = beta(i,A) + L_A*(t(1)-t_prev)*w_A';    % Compute coefficient estimates corresponding to a specific t
            t(1) = [];
            ii = ii+1;
        end
        t_prev = t_now;
    end
    
    mu = mu_old + gamma(i)*u_A; % Update mean vector
    mu_old = mu;
    beta = [beta; beta_tmp'];
end

% If smallest condition boundary-value defined & reached,
% output only the results matching these condition value(s)
if ii > 1,
    noCons = (tt > norm(beta_tmp,1));
    if 0 < sum(noCons),
        beta_t(noCons,:) = repmat(beta_tmp',sum(noCons),1);
    end
    beta = beta_t;
end

%% ____________ %% FUNCTIONS %% ___________ %%

% Normalize columns of X to have mean zero and length one.
function sX = normalize(X)

[n,p] = size(X);
sX = X-repmat(mean(X),n,1);
sX = sX*diag(1./sqrt(ones(1,n)*sX.^2));


% Find the minimum and its index over the (strictly) positive part of X
% matrix
function [m, I, J] = minplus(X)

% Remove complex elements and reset to Inf
[I,J] = find(0~=imag(X));
for i = 1:length(I),
    X(I(i),J(i)) = Inf;
end

X(X<=0) = Inf;
m = min(min(X));
[I,J] = find(X==m);
