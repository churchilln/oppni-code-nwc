function nii = load_untouch_nii(path)

    nii = struct();
    vol = niftiread(path);
    info = niftiinfo(path);

    nii.img = vol;
    nii.info = info;

    nii.hdr = struct();
    nii.hdr.dime = info.raw;
    nii.hdr.hist = info.raw;

return