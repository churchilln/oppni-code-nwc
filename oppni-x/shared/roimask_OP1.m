function roimask_OP1( pcsf_map, pwm_map, pgm_map, tsd_map, brain_mask, odir, ParamCell )
%
% functional roi masking..
%
% odir => [outpath,'/group_level/masks']
%
% --> checks for consistency of roimask!

if nargin<7 || isempty(ParamCell) || (numel(ParamCell)==1 && isempty(ParamCell{1}))
    spec_set=struct();
    spec_set.gm_smo      = 0.85;
    spec_set.csf_pthresh = 0.90;
    spec_set.wm_pthresh  = 0.90;
    spec_set.gm_rthresh  = 0.33;
    spec_set.tsd_rthresh = 0.90;
else
    spec_set.gm_smo      = str2num( ParamCell{1} );
    spec_set.csf_pthresh = str2num( ParamCell{2} );
    spec_set.wm_pthresh  = str2num( ParamCell{3} );
    spec_set.gm_rthresh  = str2num( ParamCell{4} );
    spec_set.tsd_rthresh = str2num( ParamCell{5} );
end

Mw = load_untouch_niiz(brain_mask);
mask_brain = double(Mw.img);

% masking the masks...
%
% CSF -- conservative mask, to avoidd vascular/GM spillover; no smooth since it can wipe out probs 
Vw=load_untouch_niiz(pcsf_map);
tmpvol = smooth3( Vw.img, 'gaussian',[7 7 7], spec_set.gm_smo);
volvec(:,1) = tmpvol(mask_brain>0);
Vw.img = clust_up( double( Vw.img > spec_set.csf_pthresh ) ,20);

if exist([odir,'/func_CSF_mask_grp.nii'],'file')
    Vold = load_untouch_niiz([odir,'/func_CSF_mask_grp.nii']);
    if sum(abs(Vold.img(:)-Vw.img(:)))>eps
        error('roimask - csf mask seems to have changed  - either restore param file, or delete P2 /groupmasks and try again with new setting!');
    end
else
    save_untouch_niiz(Vw,[odir,'/func_CSF_mask_grp.nii']);
end
clear Vw;


% WM -- conservative to avoid GM; smooth to wipe out peripheral WM with lots of GM pve  
Vw=load_untouch_niiz(pwm_map);
tmpvol = smooth3( Vw.img, 'gaussian',[7 7 7], spec_set.gm_smo);
volvec(:,2) = tmpvol(mask_brain>0);
Vw.img = smooth3( double(Vw.img), 'gaussian',[7 7 7], spec_set.gm_smo);
Vw.img = clust_up( double( Vw.img > spec_set.wm_pthresh ) ,20);

if exist([odir,'/func_WM_mask_grp.nii'],'file')
    Vold = load_untouch_niiz([odir,'/func_WM_mask_grp.nii']);
    if sum(abs(Vold.img(:)-Vw.img(:)))>eps
        error('roimask - wm mask seems to have changed  - either restore param file, or delete P2 /groupmasks and try again with new setting!');
    end
else
    save_untouch_niiz(Vw,[odir,'/func_WM_mask_grp.nii']);
end
clear Vw;


% GM
Vw=load_untouch_niiz(pgm_map);
tmpvol = smooth3( Vw.img, 'gaussian',[7 7 7], spec_set.gm_smo);
volvec(:,3) = tmpvol(mask_brain>0);
tmp = mask_brain; tmp(tmp>0)= double( (volvec(:,3) ./ sum(volvec,2)) >= spec_set.gm_rthresh );
Vw.img = tmp; clear tmp;

if exist([odir,'/func_GM_mask_grp.nii'],'file')
    Vold = load_untouch_niiz([odir,'/func_GM_mask_grp.nii']);
    if sum(abs(Vold.img(:)-Vw.img(:)))>eps
        error('roimask - gm mask seems to have changed  - either restore param file, or delete P2 /groupmasks and try again with new setting!');
    end
else
    save_untouch_niiz(Vw,[odir,'/func_GM_mask_grp.nii']);
end
clear Vw;


% Var
Vw=load_untouch_niiz(tsd_map);
Vw.img = clust_up( double( Vw.img > prctile(Vw.img(mask_brain>0), round(100*spec_set.tsd_rthresh) ) ) .* mask_brain ,20);
save_untouch_niiz(Vw,[odir,'/func_tSD_mask_grp.nii']);

if exist([odir,'/func_tSD_mask_grp.nii'],'file')
    Vold = load_untouch_niiz([odir,'/func_tSD_mask_grp.nii']);
    if sum(abs(Vold.img(:)-Vw.img(:)))>eps
        error('roimask - tsd mask seems to have changed  - either restore param file, or delete P2 /groupmasks and try again with new setting!');
    end
else
    save_untouch_niiz(Vw,[odir,'/func_tSD_mask_grp.nii']);
end
clear Vw;





%     % masking the masks...
%     %
%     % CSF -- conservative mask, to avoidd vascular/GM spillover; no smooth since it can wipe out probs 
%     Vw=load_untouch_niiz([outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_pCSF_grp.nii']);
%     tmpvol = smooth3( Vw.img, 'gaussian',[7 7 7], 0.85);
%     volvec(:,1) = tmpvol(mask_brain>0);
%     Vw.img = clust_up( double( Vw.img > 0.90 ) ,20);
%     mask_csf = double(Vw.img);
%     save_untouch_niiz(Vw,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_CSF_mask_grp.nii']);
%     % WM -- conservative to avoid GM; smooth to wipe out peripheral WM with lots of GM pve  
%     Vw=load_untouch_niiz([outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_pWM_grp.nii']);
%     tmpvol = smooth3( Vw.img, 'gaussian',[7 7 7], 0.85);
%     volvec(:,2) = tmpvol(mask_brain>0);
%     Vw.img = smooth3( double(Vw.img), 'gaussian',[7 7 7], 0.85);
%     Vw.img = clust_up( double( Vw.img > 0.90 ) ,20);
%     mask_wm = double(Vw.img);
%     save_untouch_niiz(Vw,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_WM_mask_grp.nii']);
%     % GM
%     Vw=load_untouch_niiz([outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_pGM_grp.nii']);
%     tmpvol = smooth3( Vw.img, 'gaussian',[7 7 7], 0.85);
%     volvec(:,3) = tmpvol(mask_brain>0);
%     tmp = mask_brain; tmp(tmp>0)= double( (volvec(:,3) ./ sum(volvec,2)) >=0.33 );
%     Vw.img = tmp; clear tmp;
%     mask_gm = double(Vw.img);
%     save_untouch_niiz(Vw,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_GM_mask_grp.nii']);
%     % Var
%     Vw=load_untouch_niiz([outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_tSD_grp.nii']);
%     Vw.img = clust_up( double( Vw.img > prctile(Vw.img(mask_brain>0),90) ) .* mask_brain ,20);
%     mask_tsd = double(Vw.img);
%     save_untouch_niiz(Vw,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_tSD_mask_grp.nii']);