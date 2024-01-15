function out = hurst_wavelet( X, j1, j2 )
%
% Wavelet monofractal estimation -- uses method of Aubry et al.
%
% syntax:
%           Hurst = hurst_wavelet(X)
%
% inputs:
%           X = data matrix, arranged (samples x time)
%
% outputs:
%           Hurst = vector of hurst exponent values


%[alphaest,cfCest,cfest,Cest,Q,j1opt,yj,varj] = LDestimate(data,regu,j1,j2,discrete_init,calcj1,printout)
%
%  Input:   data:  the input data as a row vector:  this is in fact the sequence of wavelet "approximation"
%                  coefficients, or an approximation thereof (for example often sampled data is used here).
%           regu:  (regularity) =  number of vanishing moments (of the Daubechies wavelet).
%           j1:    the lower limit of the scales chosen,  1<= j1 <= scalemax-1
%           j2:    the upper limit of the octaves chosen, 2<= j2 <= scalemax
%           discrete_init:   1:  perform the special MRA initialisation for intrinsically discrete series
%                         else: assume input data already initialised (ie is already the approximation sequence)
%           calcj1:  1: decision to run the newchoosej1 function, which will output a plot of Q(j_1) 
%                    vs j_1, as well as returning the optimal j1 values (one per j2 input to it).
%                    The values of j2 chosen are not input parameters- need to edit this function.
%                    NOTE: choosing this does Not affect the j1 passed to
%                    LDestimate for estimation. 
%
regu = 3;
% j1   = 2;
% j2   = 4;


for(i=1:size(X,1))
%     disp(['running ',num2str(i),'/',num2str(size(X,1))]);
    [Alpha(i,1),cfC,cf,C,Q,j1opt,yj(i,:),varj] = LDestimate( X(i,:)' ,regu,j1,j2,1,1,0);
end
% convert scaling exponent to Hurst exponent
out.Hurst = (Alpha+1)./2;
out.yj    = yj(:,j1:j2);
out.scal  = j1:j2;

% figure, boxplot (Hurst);
% ylabel('hurst exponent values'); title('monofractal wavelet scaling'); ylim([0 1]);
% 
