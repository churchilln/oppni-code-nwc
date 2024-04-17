function awarp_AN1(Adataset, Maskdset, Basedset, odir, ParamCell)

% Adataset = 'inputdataset.nii';  %# req/ input dataset
% SubID    = 'subjID';            %# req/ the subject ID
% Basedset = 'basedset.nii';      %# req/ reference dset- must have 4 bricks
% odir     = 'output';            %# opt/ output dir

% based on: https://github.com/ANTsX/ANTs/blob/master/Scripts/newAntsExample.sh
%
% error('This option currently disabled -- still needs detailed QC!')


pref = [odir,'/__opptmp_p2anat_warp'];

if ~exist( sprintf('%s/anat_warped.nii.gz',odir) ,'file')
    
    unix(sprintf('mkdir -p %s',pref));

    % check for path  to base file, create tpath
    if exist(Basedset,'file')
        [tpath,~,~] = fileparts(Basedset);
    else
        error('cannot find path to Based file');
    end
    
    % Require it to have enough bricks to really be a ref
    [~,nvolbase] = unix(sprintf('3dinfo -nv %s',Basedset)); nvolbase = str2num(nvolbase);
    if nvolbase < 4
        error('needs 4 volumes to serve as base!?!?!');
    end


    % make a copy
    if contains(Adataset,'.nii.gz')
        unix(sprintf('cp %s %s/anat.nii.gz',Adataset,pref));
    elseif contains(Adataset,'.nii')
        unix(sprintf('cp %s %s/anat.nii',Adataset,pref));
        unix(sprintf('gzip %s/anat.nii',pref));
    else
        error('unregocnized Adataset format');
    end

    % make a copy
    [~,refp, refe]=fileparts( Basedset );
    if strcmp(refe,'.nii')
        unix(sprintf('cp %s %s/Base.nii',Basedset,pref))
        unix(sprintf('gzip %s/Base.nii',pref));
    elseif strcmp(refe,'.gz')
        unix(sprintf('cp %s %s/Base.nii.gz',Basedset,pref))
    else
        error('template has unrecognized file ending!')
    end
    unix(sprintf('3dTcat -prefix %s/Base_skulloff.nii.gz ''%s/Base.nii.gz[0]''',pref,pref));
    unix(sprintf('3dTcat -prefix %s/Base_skullon.nii.gz ''%s/Base.nii.gz[1]''',pref,pref));
    unix(sprintf('3dTcat -prefix %s/Base_mask.nii.gz ''%s/Base.nii.gz[3]''',pref,pref));

    unix(sprintf('N4BiasFieldCorrection -d 3 -v 1 -s 4 -b [ 180 ] -c [ 50x50x50x50, 0.0 ] -i %s/anat.nii.gz -x %s -o [ %s/anat_N4.nii.gz, %s/anat_BiasField.nii.gz ]',pref,Maskdset,pref,pref))

    unix(sprintf('3dcalc -prefix %s/anat_N4_masked.nii.gz -a %s/anat_N4.nii.gz -b %s -expr ''a*b''',pref,pref,Maskdset));

    % params from "newAntsExample.sh"+
    fixd = sprintf('%s/Base_skulloff.nii.gz',pref);% 'MNI152_2009_template_SSW_brain.nii';
    movn = sprintf('%s/anat_N4_masked.nii.gz',pref);% sprintf('%s/amskBrainExtractionBrain.nii.gz',opath2a);
    nomen = sprintf('%s/aln',pref);
    unix(sprintf(['antsRegistration ', ...
    '-d 3 ', ...
    '-r [ %s, %s , 1 ] ', ...
    '-m mattes[  %s, %s , 1 , 32, regular, 0.3 ] ', ...
    '-t translation[ 0.1 ] ', ...
    '-c [ 10000x111110x11110,1.e-8,20 ] ', ...
    '-s 4x2x1vox  ', ...
    '-f 6x4x2 -l 1 ', ...
    '-m mattes[  %s, %s , 1 , 32, regular, 0.3 ] ', ...
    '-t rigid[ 0.1 ] ', ...
    '-c [ 10000x111110x11110,1.e-8,20 ] ', ...
    '-s 4x2x1vox  ', ...
    '-f 3x2x1 -l 1 ', ...
    '-m mattes[  %s, %s , 1 , 32, regular, 0.3 ] ', ...
    '-t affine[ 0.1 ] ', ...
    '-c [ 10000x111110x11110,1.e-8,20 ] ', ...
    '-s 4x2x1vox  ', ...
    '-f 3x2x1 -l 1 ', ...
    '-m mattes[  %s, %s , 0.5 , 32 ] ', ...
    '-m cc[  %s, %s , 0.5 , 4 ] ', ...
    '-t SyN[ .20, 3, 0 ] ', ...
    '-c [ 100x100x50,-0.01,5 ] ', ...
    '-s 1x0.5x0vox ', ...
    '-f 4x2x1 -l 1 -u 1 -z 1 ', ...
    '-o [ %s,%s_diff.nii.gz,%s_inv.nii.gz]'],...
    fixd, movn, fixd, movn, fixd, movn, fixd, movn, fixd, movn, fixd, movn, nomen, nomen, nomen ));
    %
    unix(sprintf('antsApplyTransforms -d 3 -i %s -r %s -n linear -t %s1Warp.nii.gz -t %s0GenericAffine.mat -o %s_warped.nii.gz',movn,fixd,nomen,nomen,nomen));

    %--transfer over the relevant files...
    unix(sprintf('cp %s/anat_N4.nii.gz %s/anat_proc.nii.gz',pref,odir));
    unix(sprintf('cp %s/anat_N4_masked.nii.gz %s/anat_procss.nii.gz',pref,odir));
    unix(sprintf('cp %s/aln_warped.nii.gz %s/anat_warped.nii.gz',pref,odir));
    unix(sprintf('cp %s/aln1Warp.nii.gz %s/anatQQ_WARP.nii.gz',pref,odir));
    unix(sprintf('cp %s/aln0GenericAffine.mat %s/anatQQ_GenericAffine.mat',pref,odir));
    
    unix(sprintf('rm -rf %s',pref));
else
    disp('ants-warp already exists!')
end