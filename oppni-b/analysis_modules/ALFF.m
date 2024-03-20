function output = ALFF( datacell, params )
%
% =========================================================================
% MODULE_FALFF: fractional amplitude of low-frequency fluctuations
% =========================================================================
%
%   Syntax:
%           output = module_alff( datamat, split_info )
%
% ------------------------------------------------------------------------%

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'ALFF';
output.attributes.design_type   = 'nocontrast';
output.attributes.model_type    = 'univariate';
output.attributes.image_types   = 'ALFF,fALFF,zALFF,zfALFF';
output.attributes.mat2d_types   = [];
output.attributes.stats_types   = [];
output.attributes.uses_roifile  = 0;
output.attributes.uses_taskfile = 0;

if(nargin==0)
    disp('no inputs - returning attributes');
    return;
end
%----------------------- Default Parameter Checks ------------------------%
%-------------------------------------------------------------------------%

Nr = numel(datacell);

for ir=1:Nr

    datamat = datacell{ir};

    % dimensions
    [Nvox Ntime] =  size(datamat); 
    TR           = round(params.TR_MSEC./1000);
    % parameters for spectral estimation --> splithalf not yet used
    %Nhalf        = floor(Ntime/2);
    %Fny    = 0.5 * (1/TR);                % nyquist frequency
    %NFFT   = 2^nextpow2( Nhalf );         % Next power of 2 from length of time-axis
    %f      = Fny*linspace(0,1,NFFT/2+1);  % fourier data corresponds to these frequency points
    
    % welch power estimation
    [powSum,f] = pwelch_matrix( datamat,[],1/TR,['hann'],[ 2^(nextpow2(Ntime)-1) ],[] );
    ampSum = sqrt(powSum); clear powSum; %- convert power to "amplitude" spectrum
    
    % frequency maps
    output.image.ALFF(:,ir)   = mean(ampSum(:,f>=0.01 & f<=0.08),2);
    output.image.fALFF(:,ir)  =  sum(ampSum(:,f>=0.01 & f<=0.08),2) ./ (sum(ampSum(:,f>0),2)+eps);
    output.image.zALFF(:,ir)  = zscore( output.image.ALFF(:,ir)  );
    output.image.zfALFF(:,ir) = zscore( output.image.fALFF(:,ir) );
end

% averaging across runs
output.image.ALFF(:,ir)  = mean(output.image.ALFF(:,ir),3);
output.image.fALFF(:,ir) = mean(output.image.fALFF(:,ir),3);
output.image.zALFF(:,ir) = mean(output.image.zALFF(:,ir),3);
output.image.zfALFF(:,ir) = mean(output.image.zfALFF(:,ir),3);
