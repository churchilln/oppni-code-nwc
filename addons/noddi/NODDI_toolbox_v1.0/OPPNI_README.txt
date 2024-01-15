modified components of code responsible for NIFTI I/O, to avoid dependency
on SPM toolbox (niftimatlab), to ensure (a) consistency in data formatting,
and (b) to avoid issues with MEX functions.

modified code is listed below. all changes are marked in the scripts with
an [NWC ed] tag.

CreateROI.m
SaveAsNIfTI.m
SaveParamsAsNIfTI.m

