%
% -------------------------------------------------------------------------
%          module: ROIREG
% -------------------------------------------------------------------------
%  pipeline stage: P2
%      applied to: fMRI data
%         summary: regressed out signal from WM/CSF noise ROIs, correction for physio noise 
%   function call: [Xreg,stat] = roireg_<name>( functional_run, roi_paths, ParamCell )
%                 functional_run = path+name of functional data signal is being extacted from 
%                      roi_paths = cell array containing paths to {1} brain masks, {2} parcellations 
%                     (ParamCell = supplemental parameters)
%
