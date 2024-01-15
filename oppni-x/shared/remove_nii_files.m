function remove_nii_files( infile, outfile, idx_tokeep, strip_hd )
%
% for "scrubbing", i.e. deleting specific files from a 3D functional run
%

if nargin<4
    strip_hd=0;
else
    if contains(outfile,'.nii.gz')
        outfile_unfix = [outfile(1:end-7),'_unfix.nii.gz'];
    elseif contains(outfile,'.nii')
        outfile_unfix = [outfile(1:end-4),'_unfix.nii'];
    else
        error('unknown infile format')
    end
end

Vi = load_untouch_niiz(infile);
Nt = size(Vi.img,4);

TMPVOL = Vi.img;
TMPVOL(:,:,:,idx_tokeep==0)=[]; % binary label

Vo=Vi;
Vo.img = TMPVOL;
Vo.hdr.hist = Vi.hdr.hist;
Vo.hdr.dime.dim(5) = size(TMPVOL,4);

if strip_hd==0
    save_untouch_niiz(Vo,outfile);
else
    save_untouch_niiz(Vo,outfile_unfix);
    unix(sprintf('nifti_tool -rm_ext ALL -prefix %s -infiles %s',outfile,outfile_unfix));
    if exist(outfile,'file')
        unix(sprintf('rm %s',outfile_unfix))
    else
        error('failed to create new file w stripped header')
    end
end
