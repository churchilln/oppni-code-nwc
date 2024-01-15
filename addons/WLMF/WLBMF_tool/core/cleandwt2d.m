function [coefs,jssbord]=cleandwt2d(data,Nwt,norm,jmax);
% function [coefs,jssbord]=cleandwt2d(data,Nwt,norm,jmax);
%
% Stephane Roux, Lyon, 2004
% Herwig Wendt, Lyon, 2006 - 2008 (last modification)

J = log2(length(data)) ;
if nargin < 3; jmax=J; end
    
decompnorm=-.5;
expo=norm-decompnorm;

h = daubcqf(2*Nwt,'min');

Lh=length(h);
wt = abs(mdwt(data,h,jmax));

for j=1:jmax
    % we normalize the coefs
    coefs(j).x  = [wt( 1:2^(J-j)             , (2^(J-j) +1):2^(J-j+1))] * 2^(j*2*expo);
    coefs(j).y  = [wt((2^(J-j) +1):2^(J-j+1) , 1:2^(J-j))] * 2^(j*2*expo);
    coefs(j).xy = [wt((2^(J-j) +1):2^(J-j+1) , (2^(J-j) +1):2^(J-j+1))] * 2^(j*2*expo);
    
    % take away border coefficients
    Nbredund(j)=(Lh-1)*(2^j-1);        % borders of redundant wavelet transform
    Nbord(j)=floor(Nbredund(j)/2^j);   % borders of dyadic true
    if length(coefs(j).x) > Nbord(j)
        jssbord=j-1;
        coefs(j).xssb=coefs(j).x(1:end-Nbord(j),1:end-Nbord(j));
        coefs(j).yssb=coefs(j).y(1:end-Nbord(j),1:end-Nbord(j));
        coefs(j).xyssb=coefs(j).xy(1:end-Nbord(j),1:end-Nbord(j));
    end
end
end