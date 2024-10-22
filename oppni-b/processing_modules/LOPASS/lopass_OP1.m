function datamat_filt = lopass_OP1( datamat, TR_MSEC, ParamCell )
%
% .lopass_OP1:
% .low-pass filtering to remove high-frequency signal fluctuations  (mainly thermal/physiological noise) 
% .uses a Butterworth filter with >0.10 stopband / <0.08 passband

[Nvox Ntime] = size(datamat); % matrix dimensions
% using simple Butterworth filter -- linear phase/ flat frequency, rolloff not great but this is tolerable for fmri
Wp = (2*(TR_MSEC/1000))*0.08; % passband is below 0.08 Hz
Ws = (2*(TR_MSEC/1000))*0.10; % stopband is above 0.10 Hz
% filter design: max passband attn. =50% / min stopband attn =1%
[Nord, Wcut] = buttord( Wp, Ws, 3,10 );
% lowpass butterworth filter with desired cutoff
[B1,A1] = butter(Nord,Wcut,'low');
% zero-phase forward/reverse filter
datamat_filt  = filtfilt( B1,A1, datamat' )';
