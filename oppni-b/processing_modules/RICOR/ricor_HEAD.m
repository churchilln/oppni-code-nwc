%
% -------------------------------------------------------------------------
%          module: RICOR
% -------------------------------------------------------------------------
%  pipeline stage: P1
%      applied to: fMRI data
%         summary: regresses out physiologic noise based on cardiac & respiratory data 
%   function call: ricor_<name>( Funcfile, prefix, odir, pulsfile, respfile, acqpar, ParamCell )
%                       Funcfile = path+full name of uncorrected fMRI file
%                         prefix = prefix given to the corrected file (output)
%                           odir = directory to store corrected file (output)
%                       pulsfile = cardiac pulsation file
%                       respfile = respiration file
%                         acqpar = structure containing acquisition params (tpatt, ndrop, tr_msec, physamp_msec) 
%                     (ParamCell = supplemental parameters)
%
