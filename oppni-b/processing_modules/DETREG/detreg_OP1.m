function [Xreg,stat] = detreg_OP1( Nt, TR_MSEC, ParamCell )
%
% .detreg_OP1:
% .estimation of low-frequency basis for general denoising 
% .uses a legendre polynomial basis

ord = ParamCell{1};

if nargin<3
    TR_MSEC = NaN;
end

if ischar(ord)
    if(strcmpi(ord,'A') || strcmpi(ord,'AUTO'))
        if ~isfinite(TR_MSEC)
            error('need to specify TR_MSEC field for auto-detrending!')
        end
        ord = 1+floor( (TR_MSEC/1000) * Nt ./ 150 );
    else
        ord = str2num(ord);
        if isempty(ord)
            error('detreg input must be numeric or "A" for auto-estimation')
        end
    end
end

stat = [ord];

x  = linspace(-1,1,Nt)';
%
DL = [];
%
if (ord>=0)
    DL = [];
    for i = 0:ord
        X = generate_legendre(i,Nt);
        DL = [DL X];
    end
end

Xreg = DL;

function X = generate_legendre(ord,N)

x  = linspace(-1,1,N)';
pk = LegendrePoly(ord);
xp = bsxfun(@power,x,[ord:-1:0]);    
X = bsxfun(@times,xp,pk');    
X = sum(X,2);

function pk = LegendrePoly(n)
% for nonnegative int n, compute Legendre polynomial P_ns. 
% Return the result as a vector whose mth element is the coefficient of x^(n+1-m).
% polyval(LegendrePoly(n),x) evaluates P_n(x).

if n==0 
    pk = 1;
elseif n==1
    pk = [1 0]';
else
    pkm2 = zeros(n+1,1);
    pkm2(n+1) = 1;
    pkm1 = zeros(n+1,1);
    pkm1(n) = 1;

    for k=2:n
        
        pk = zeros(n+1,1);

        for e=n-k+1:2:n
            pk(e) = (2*k-1)*pkm1(e+1) + (1-k)*pkm2(e);
        end
        
        pk(n+1) = pk(n+1) + (1-k)*pkm2(n+1);
        pk = pk/k;
        
        if k<n
            pkm2 = pkm1;
            pkm1 = pk;
        end
    end
end

% if( ord>=0 ) d0 = ones(N,1);                                       DL=[DL d0]; end
% if( ord>=1 ) d1 = x;                                               DL=[DL d1]; end
% if( ord>=2 ) d2 = 0.5    *(3*x.^2   - 1);                          DL=[DL d2]; end
% if( ord>=3 ) d3 = 0.5    *(5*x.^3   - 3*x);                        DL=[DL d3]; end
% if( ord>=4 ) d4 = 0.125  *(35*x.^4  - 30*x.^2  + 3);               DL=[DL d4]; end
% if( ord>=5 ) d5 = 0.125  *(63*x.^5  - 70*x.^3  + 15*x);            DL=[DL d5]; end
% if( ord>=6 ) d6 = 0.0625 *(231*x.^6 - 315*x.^4 + 105*x.^2 - 5);    DL=[DL d6]; end
% if( ord>=7 ) d7 = 0.0625 *(429*x.^7 - 693*x.^5 + 315*x.^3 - 35*x); DL=[DL d7]; end
