%
% -------------------------------------------------------------------------
%          module: MOTREG
% -------------------------------------------------------------------------
%  pipeline stage: P2
%      applied to: fMRI data
%         summary: regressed out signal from motion parameter estimates, correction for motion noise 
%   function call: [Xreg,stat] = motreg_<name>( mpefile, ParamCell )
%                      mpefile = location + name of motion parameter estimate file
%                     (ParamCell = supplemental parameters)
%
