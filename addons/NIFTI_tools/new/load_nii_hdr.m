function hdr = load_nii_hdr(path)

    %nii = struct();
    %vol = niftiread(path);
    info = niftiinfo(path);

    %nii.img = vol;
    %nii.info = info;

    hdr = struct();
    hdr.dime = info.raw;
    hdr.hist = info.raw;

return