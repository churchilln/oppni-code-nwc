function [acav,acd1,acd2] = estim_acd_temp( img, mask )
%
%   temporal autocorrelation estimator
%
%   [aca,acd1,acd2] = estim_acd_temp( img, mask )
%
%   where aca  = autocorrelation average
%         acd1 = autocorrelation difference, 1st order
%         acd2 = autocorrelation difference, 2nd order
%

for lag=1:3
    % taxis
    a= img(:,:,:,1:end-lag); 
    b= img(:,:,:,lag+1:end);
    c= mask;
    %
    a= (a-mean(a,4))./(std(a,0,4)+eps);
    b= (b-mean(b,4))./(std(b,0,4)+eps);
    tmp = sum(a.*b.*c,4)./(size(a,4)-1);
    stt(lag,1) = mean(tmp(c>0));
end
acav = mean(stt);
acd1 = mean( stt(2:end) - stt(1:end-1) );
acd2 = mean( stt(3:end) - 2*stt(2:end-1) + stt(1:end-2) );

% technically, acd should be del / TR-size to give it scale in mm
