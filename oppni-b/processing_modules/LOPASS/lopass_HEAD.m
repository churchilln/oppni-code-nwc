%
% -------------------------------------------------------------------------
%          module: LOPASS
% -------------------------------------------------------------------------
%  pipeline stage: P2
%      applied to: fMRI data
%         summary: filtering of high-frequency physio signal - correction for physio noise 
%   function call: datamat_filt = roireg_<name>( datamat, TR_MSEC, ParamCell )
%                        datamat = matricized functional data, as a voxel x time matrix 
%                        TR_MSEC = repetition time, in milliseconds 
%                     (ParamCell = supplemental parameters)
%
