function save_untouch_niiz(nii, path, mni_space)
    if ~exist("mni_space","var")
        mni_space = false;
    end

    if contains(path, ".nii.gz")
        save_untouch_nii(nii, path, true, mni_space);
    else
        save_untouch_nii(nii, path, false, mni_space);
    end
return