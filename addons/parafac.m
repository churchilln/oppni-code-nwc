function [A,B,C,fit,it] = parafac(X,DimX,Fac,crit,Constraints,A,B,C,maxit,DoLineSearch);

% Complex PARAFAC-ALS
% Fits the PARAFAC model Xk = A*Dk*B.' + E
% where Dk is a diagonal matrix holding the k'th
% row of C.
%
% Uses on-the-fly projection-compression to speed up 
% the computations. This requires that the first mode 
% is the largest to be effective
% 
% INPUT
% X          : Data
% DimX       : Dimension of X
% Fac        : Number of factors
% OPTIONAL INPUT
% crit       : Convergence criterion (default 1e-6)
% Constraints: [a b c], if e.g. a=0 => A unconstrained, a=1 => A nonnegative
% A,B,C      : Initial parameter values
%
% I/O
% [A,B,C,fit,it]=parafac(X,DimX,Fac,crit,A,B,C);
%
% Copyright 1998
% Rasmus Bro
% KVL, Denmark, rb@kvl.dk

% Initialization
if nargin<9
  maxit   = 2500;      % Maximal number of iterations
end
showfit = pi;         % Show fit every 'showfit'th iteration (set to pi to avoid)

if nargin<4
  crit=1e-6;
end

if crit==0
  crit=1e-6;
end

I = DimX(1);
J = DimX(2);
K = DimX(3);

InitWithRandom=0;
if nargin<8
   InitWithRandom=1;
end
if nargin>7 & size(A,1)~=I
  InitWithRandom=1;
end

if nargin<5
   ConstA = 0;ConstB = 0;ConstC = 0;
else
   ConstA = Constraints(1);ConstB = Constraints(2);ConstC = Constraints(3);
end

if InitWithRandom

  if I<Fac
    A = rand(I,Fac);
  else
    A = orth(rand(I,Fac));
  end
  if J<Fac
    B = rand(J,Fac);
  else
    B = orth(rand(J,Fac));
  end
  if K<Fac
    C = rand(K,Fac);
  else
    C = orth(rand(K,Fac));
  end
end

SumSqX = sum(sum(abs(X).^2));
fit    = SumSqX;
fit0   = fit;
fitold = 2*fit;
it     = 0;
Delta  = 5;

while abs((fit-fitold)/fitold)>crit&it<maxit&fit>10*eps
   it=it+1;
   fitold=fit;

   % Do line-search
   if rem(it+2,2)==-1
      [A,B,C,Delta]=linesrch(X,DimX,A,B,C,Ao,Bo,Co,Delta);
   end
   
   Ao=A;Bo=B;Co=C;
   % Update A
   Xbc=0;
   for k=1:K
     Xbc = Xbc + X(:,(k-1)*J+1:k*J)*conj(B*diag(C(k,:)));
   end
   if ConstA == 0 % Unconstrained
      A = Xbc*pinv((B'*B).*(C'*C)).';
   elseif ConstA == 1 % Nonnegativity, requires reals
      Aold = A;
      for i = 1:I
         ztz = (B'*B).*(C'*C);
         A(i,:) = fastnnls(ztz,Xbc(i,:)')';
      end
      if any(sum(A)<100*eps*I)
         A = .99*Aold+.01*A; % To prevent a matrix with zero columns
      end
   elseif ConstA == 2 % Orthogonality
      A = Xbc*(Xbc'*Xbc)^(-.5);
   elseif ConstA == 3 % Unimodality
      A = unimodalcrossproducts((B'*B).*(C'*C),Xbc',A);
   end

   % Project X down on orth(A) - saves time if first mode is large
   [Qa,Ra]=qr(A,0);
   x=Qa'*X;

   % Update B
   if ConstB == 10 % Procrustes
      B = eye(Fac);
   else
      Xac=0;
      for k=1:K
         Xac = Xac + x(:,(k-1)*J+1:k*J).'*conj(Ra*diag(C(k,:)));
      end
      if ConstB == 0 % Unconstrained
         B = Xac*pinv((Ra'*Ra).*(C'*C)).';
      elseif ConstB == 1 % Nonnegativity, requires reals
         Bold = B;
         for j = 1:J
            ztz = (Ra'*Ra).*(C'*C);
            B(j,:) = fastnnls(ztz,Xac(j,:)')';
         end
         if any(sum(B)<100*eps*J)
            B = .99*Bold+.01*B; % To prevent a matrix with zero columns
         end
      end
   end
  
    % Update C
    if ConstC == 0 % Unconstrained
       ab=pinv((Ra'*Ra).*(B'*B));
       for k=1:K 
          C(k,:) = (ab*diag(Ra'* x(:,(k-1)*J+1:k*J)*conj(B))).';
       end
    elseif ConstC == 1  % Nonnegativity, requires reals
       Cold = C;
       ztz = (Ra'*Ra).*(B'*B);
       for k = 1:K
          xab = diag(Ra'* x(:,(k-1)*J+1:k*J)*B);
          C(k,:) = fastnnls(ztz,xab)';
       end
       if any(sum(C)<100*eps*K)
          C = .99*Cold+.01*C; % To prevent a matrix with zero columns
       end
    elseif ConstC == 2 % Orthogonality
       Z=(Ra'*Ra).*(B'*B);
       Y=[];
       for k=1:K
          d=diag(Ra'*x(:,(k-1)*J+1:k*J)*B)'; 
          Y=[Y;d];
       end;
       [P,D,Q]=svd(Y,0);
       C=P*Q';
    elseif ConstC == 3 % Unimodality
       xab = [];
       for k = 1:K
          xab = [xab diag(Ra'* x(:,(k-1)*J+1:k*J)*B)];
       end
       C = unimodalcrossproducts((Ra'*Ra).*(B'*B),xab,C);
    elseif ConstC == 10 % GPA => Isotropic scaling factor
       ab=(Ra'*Ra).*(B'*B);
       ab = pinv(ab(:));
       C(1,:) = 1;
       for k=2:K 
          yy = [];
          yyy = diag(Ra'* x(:,(k-1)*J+1:k*J)*conj(B)).';
          for f=1:Fac
             yy = [yy;yyy(:)];
          end
          C(k,:) = ab*yy;
       end
    end
      
    % Calculating fit. Using orthogonalization instead
   %fit=0;for k=1:K,residual=X(:,(k-1)*J+1:k*J)-A*diag(C(k,:))*B.';fit=fit+sum(sum((abs(residual).^2)));end
   [Qb,Rb]=qr(B,0);
   [Z,Rc]=qr(C,0);
   fit=SumSqX-sum(sum(abs(Ra*ppp(Rb,Rc).').^2));
   
   if rem(it,showfit)==0
      fprintf(' %12.10f       %g        %3.4f \n',fit,it,100*(1-fit/fit0));
   end
end

% ORDER ACCORDING TO VARIANCE
Tuck     = diag((A'*A).*(B'*B).*(C'*C));
[out,ID] = sort(Tuck);
A        = A(:,ID);
if ConstB ~= 10 % Else B is eye
   B        = B(:,ID);
end
C        = C(:,ID);
% NORMALIZE A AND C (variance in B)
if ConstB ~= 10 % Then B is eye
   for f=1:Fac,normC(f) = norm(C(:,f));end
   for f=1:Fac,normA(f) = norm(A(:,f));end
   B        = B*diag(normC)*diag(normA);  
   A        = A*diag(normA.^(-1));
   C        = C*diag(normC.^(-1));
   
   % APPLY SIGN CONVENTION
   SignA = sign(sum(sign(A))+eps);
   SignC = sign(sum(sign(C))+eps);
   A = A*diag(SignA);
   C = C*diag(SignC);
   B = B*diag(SignA)*diag(SignC);
end


function AB=ppp(A,B);

% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
%
% Copyright, 1998 - 
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. Furthermore, the
% code can not be made part of anything but the 'N-way Toolbox'.
% In case of doubt, contact the holder of the copyrights.
%
% Rasmus Bro
% Chemometrics Group, Food Technology
% Department of Food and Dairy Science
% Royal Veterinary and Agricultutal University
% Rolighedsvej 30, DK-1958 Frederiksberg, Denmark
% Phone  +45 35283296
% Fax    +45 35283245
% E-mail rb@kvl.dk
%
% The parallel proportional profiles product - triple-P product
% For two matrices with similar column dimension the triple-P product
% is ppp(A,B) = [kron(B(:,1),A(:,1) .... kron(B(:,F),A(:,F)]
% 
% AB = ppp(A,B);
%
% Copyright 1998
% Rasmus Bro
% KVL,DK
% rb@kvl.dk

[I,F]=size(A);
[J,F1]=size(B);

if F~=F1
   error(' Error in ppp.m - The matrices must have the same number of columns')
end

AB=zeros(I*J,F);
for f=1:F
   ab=A(:,f)*B(:,f).';
   AB(:,f)=ab(:);
end

function [x,w] = fastnnls(XtX,Xty,tol)
%NNLS	Non-negative least-squares.
%	b = fastnnls(XtX,Xty) returns the vector b that solves X*b = y
%	in a least squares sense, subject to b >= 0, given the inputs
%       XtX = X'*X and Xty = X'*y.
%
%	A default tolerance of TOL = MAX(SIZE(X)) * NORM(X,1) * EPS
%	is used for deciding when elements of b are less than zero.
%	This can be overridden with b = fastnnls(X,y,TOL).
%
%	[b,w] = fastnnls(XtX,Xty) also returns dual vector w where
%	w(i) < 0 where b(i) = 0 and w(i) = 0 where b(i) > 0.
%
%	See also LSCOV, SLASH.

%	L. Shure 5-8-87
%	Revised, 12-15-88,8-31-89 LS.
%	Copyright (c) 1984-94 by The MathWorks, Inc.

%       Revised by:
%	Copyright
%	Rasmus Bro 1995
%	Denmark
%	E-mail rb@kvl.dk
%       According to Bro & de Jong, J. Chemom, 1997

% initialize variables


if nargin < 3
    tol = 10*eps*norm(XtX,1)*max(size(XtX));
end
[m,n] = size(XtX);
P = zeros(1,n);
Z = 1:n;
x = P';
ZZ=Z;
w = Xty-XtX*x;

% set up iteration criterion
iter = 0;
itmax = 30*n;

% outer loop to put variables into set to hold positive coefficients
while any(Z) & any(w(ZZ) > tol)
    [wt,t] = max(w(ZZ));
    t = ZZ(t);
    P(1,t) = t;
    Z(t) = 0;
    PP = find(P);
    ZZ = find(Z);
    nzz = size(ZZ);
    z(PP')=(Xty(PP)'/XtX(PP,PP)');
    z(ZZ) = zeros(nzz(2),nzz(1))';
    z=z(:);
% inner loop to remove elements from the positive set which no longer belong

    while any((z(PP) <= tol)) & iter < itmax

        iter = iter + 1;
        QQ = find((z <= tol) & P');
        alpha = min(x(QQ)./(x(QQ) - z(QQ)));
        x = x + alpha*(z - x);
        ij = find(abs(x) < tol & P' ~= 0);
        Z(ij)=ij';
        P(ij)=zeros(1,max(size(ij)));
        PP = find(P);
        ZZ = find(Z);
        nzz = size(ZZ);
        z(PP)=(Xty(PP)'/XtX(PP,PP)');
        z(ZZ) = zeros(nzz(2),nzz(1));
        z=z(:);
    end
    x = z;
    w = Xty-XtX*x;
end

x=x(:);

