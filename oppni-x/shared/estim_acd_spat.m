function [acav,acd1,acd2] = estim_acd_spat( img, mask )
%
%   spatial autocorrelation estimator
%
%   [aca,acd1,acd2] = estim_acd_spat( img, mask )
%
%   where aca  = autocorrelation average
%         acd1 = autocorrelation difference, 1st order
%         acd2 = autocorrelation difference, 2nd order
%

for lag=1:3
    % xaxis
    a= img(1:end-lag,:,:); 
    b= img(lag+1:end,:,:); 
    c = (mask(1:end-lag,:,:) .* mask(lag+1:end,:,:)); 
    sxyz(lag,1) = sum( zscore(double(a(c>0))).*zscore(double(b(c>0))))./(sum(c(:))-1); % mean acorr
    % yaxis
    a= img(:,1:end-lag,:); 
    b= img(:,lag+1:end,:); 
    c = (mask(:,1:end-lag,:) .* mask(:,lag+1:end,:)); 
    sxyz(lag,2) = sum( zscore(double(a(c>0))).*zscore(double(b(c>0))))./(sum(c(:))-1); % mean acorr
    % zaxis
    a= img(:,:,1:end-lag); 
    b= img(:,:,lag+1:end); 
    c = (mask(:,:,1:end-lag) .* mask(:,:,lag+1:end)); 
    sxyz(lag,3) = sum( zscore(double(a(c>0))).*zscore(double(b(c>0))))./(sum(c(:))-1); % mean acorr
end
sxyz=mean(sxyz,2);

acav = mean(sxyz);
acd1 = mean( sxyz(2:end) - sxyz(1:end-1) );
acd2 = mean( sxyz(3:end) - 2*sxyz(2:end-1) + sxyz(1:end-2) );

% technically, acd should be del / vox-size to give it scale in mm
