%
% -------------------------------------------------------------------------
%          module: DETREG
% -------------------------------------------------------------------------
%  pipeline stage: P2
%      applied to: fMRI data
%         summary: regressed out signal from motion parameter estimates, correction for motion noise 
%   function call: [Xreg,stat] = detreg_<name>( Nt, TR_MSEC, ParamCell )
%                             Nt = number of time points in scan run
%                        TR_MSEC = repetition time in milliseconds
%                     (ParamCell = supplemental parameters)
%
