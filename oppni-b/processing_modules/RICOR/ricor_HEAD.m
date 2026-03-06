%
% -------------------------------------------------------------------------
%          module: RICOR
% -------------------------------------------------------------------------
%  pipeline stage: P1
%      applied to: fMRI data
%         summary: regresses out physiologic noise based on cardiac & respiratory data 
%   function call: ricor_<name>( Funcfile, prefix, odir, physfile, acqpar, ParamCell )
%                       Funcfile = path+full name of uncorrected fMRI file
%                         prefix = prefix given to the corrected file (output)
%                           odir = directory to store corrected file (output)
%                       physfile = physiologic file string
%                         acqpar = structure containing acquisition params (tpatt, ndrop, tr_msec, physamp_msec) 
%                     (ParamCell = supplemental parameters)
%
