function save_untouch_nii(nii, path, compressed, mni_space)

    if ~exist("compressed","var")
        compressed=false;
    elseif compressed
        path = replace(path, "nii.gz", "nii");
    end

    datatypeStr = {'ubit1', 'uint8', 'int16', 'int32', 'single', 'single', ...
                        'double', 'uint8', 'int8', 'uint16', 'uint32', 'int64', ...
                        'uint64'};
    datatypeCode = {1, 2, 4, 8, 16, 16, 64, 2, 256, 512, 768, 1024, 1280};
    dataMap = containers.Map(datatypeCode, datatypeStr);
    sizeSet = {1, 8, 16, 32, 32, 32, 64, 8, 8, 16, 32, 64, 64};
    sizeMap = containers.Map(datatypeCode, sizeSet);

    nii.info.Datatype = dataMap(nii.hdr.dime.datatype);
    nii.info.BitsPerPixel = sizeMap(nii.hdr.dime.datatype);
    nii.img = cast(nii.img, nii.info.Datatype);
    niftiwrite(nii.img, path, nii.info, "Compressed",compressed);

    if exist("mni_space","var") && mni_space
        unix(sprintf('3drefit -view +tlrc -space MNI %s', path));
    end

return