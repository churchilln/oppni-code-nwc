close all;
clear;

%% This script applies spatial warping, obtained using pre-constructed AFNI+ANTs protocol, to normalize voxelwise texture maps

% hardcode absolute path
pathloc_fs = 'Neurocovid/NCOV_SUBS';
pathloc_af = 'Neurocovid/fmri_proc';
e = dir(sprintf('%s/sub-*-V01',pathloc_fs));
if isempty(e)
    error('no folders found!');
end

% %% WM WARPING
% for i=1:numel(e)
% 
%     if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
%             
%         prfix = e(i).name,
% 
%         if exist(sprintf('texturial_matfiles_vox/%s_wm_Output_Vox_restrict_texmaps_msk.nii.gz',prfix),'file')
% 
%             disp('ok!');
%         
%             warploc = sprintf('Neurocovid/fmri_proc/%s/anat_proc/subpipe_001',prfix);
%     
%             % create affine for multivolume...
%             fid = fopen(sprintf('%s/anatQQ.aff12.1D',warploc),'r'); % read in aff
%                 newline1 = fgetl(fid);
%                 newline2 = fgetl(fid);
%             fclose(fid);
%             fid = fopen(sprintf('texturial_matfiles_vox_warp_wm/%s.aff12.1D',prfix),'w'); % make new local aff
%                 fprintf(fid,'%s\n',newline1);
%                 for i=1:7
%                 fprintf(fid,'%s\n',newline2);
%                 end
%             fclose(fid);
%     
%             %%% %%% %%% %%% GM MASK FIRST
%     
%             infilea = sprintf('texturial_matfiles_vox/%s_wm_Output_Vox_restrict_texmaps_msk.nii.gz',prfix);
%             oufilea = sprintf('texturial_matfiles_vox_warp_wm/%s_wm_mask_deobed.nii.gz',prfix);
%             oufileb = sprintf('texturial_matfiles_vox_warp_wm/%s_wm_mask_reored.nii.gz',prfix);
%             oufilec = sprintf('texturial_matfiles_vox_warp_wm/%s_wm_mask_warped_2mm.nii.gz',prfix);
%     
%             % resliceded it intos the correct orientationalles
%             unix(sprintf('3dWarp -oblique2card -prefix %s -wsinc5 %s', oufilea, infilea)); %-wsinc5
%             unix(sprintf('fslreorient2std %s %s', oufilea, oufileb));        
%             % now we shall do the resamplingue of thine mapses (singlevol)
%             unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 2 -source %s -nwarp "%s/anatQQ_WARP.nii.gz %s/anatQQ.aff12.1D" -prefix %s -ainterp wsinc5',...
%                 warploc,oufileb,warploc,warploc,oufilec));
%     
%             unix(sprintf('rm %s %s',oufilea,oufileb))
%     
%             %%% %%% %%% %%% TEXTURE MAPS NEXT
%     
%             infilea = sprintf('texturial_matfiles_vox/%s_wm_Output_Vox_restrict_texmaps_av.nii.gz',prfix);
%             oufilea = sprintf('texturial_matfiles_vox_warp_wm/%s_tex_maps_deobed.nii.gz',prfix);
%             oufileb = sprintf('texturial_matfiles_vox_warp_wm/%s_tex_maps_reored.nii.gz',prfix);
%             oufilec = sprintf('texturial_matfiles_vox_warp_wm/%s_tex_maps_warped_2mm.nii.gz',prfix);
%     
%             % resliceded it intos the correct orientationalles
%             unix(sprintf('3dWarp -oblique2card -prefix %s -wsinc5 %s', oufilea, infilea)); %-wsinc5
%             unix(sprintf('fslreorient2std %s %s', oufilea, oufileb));        
%             % now we shall do the resamplingue of thine mapses (multivol)
%             unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 2 -source %s -nwarp "%s/anatQQ_WARP.nii.gz texturial_matfiles_vox_warp_wm/%s.aff12.1D" -prefix %s -ainterp wsinc5',...
%                 warploc,oufileb,warploc,prfix,oufilec));
%             
%             unix(sprintf('rm %s %s',oufilea,oufileb))
%         end
%     end
% end

% for i=1:numel(e)
%     i,
%     if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
%             
%         prfix = e(i).name,
%             oufilec = sprintf('texturial_matfiles_vox_warp_wm/%s_wm_mask_warped_2mm.nii.gz',prfix);
%             oufiled = sprintf('texturial_matfiles_vox_warp_wm/%s_wm_mask_warped_2mm_smo6.nii.gz',prfix);
%             unix(sprintf('3dmerge -prefix %s -doall -1blur_fwhm 6 %s',oufiled,oufilec));
% 
%             oufilec = sprintf('texturial_matfiles_vox_warp_wm/%s_tex_maps_warped_2mm.nii.gz',prfix);
%             oufiled = sprintf('texturial_matfiles_vox_warp_wm/%s_tex_maps_warped_2mm_smo6.nii.gz',prfix);
%             unix(sprintf('3dmerge -prefix %s -doall -1blur_fwhm 6 %s',oufiled,oufilec));
%     end
% end
% 
% %% GM WARPING
% for i=1:numel(e)
% 
%     if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
%             
%         prfix = e(i).name,
% 
%         if exist(sprintf('texturial_matfiles_vox/%s_gm_Output_Vox_restrict_texmaps_msk.nii.gz',prfix),'file')
% 
%             disp('ok!');
%         
%             warploc = sprintf('Neurocovid/fmri_proc/%s/anat_proc/subpipe_001',prfix);
%     
%             % create affine for multivolume...
%             fid = fopen(sprintf('%s/anatQQ.aff12.1D',warploc),'r'); % read in aff
%                 newline1 = fgetl(fid);
%                 newline2 = fgetl(fid);
%             fclose(fid);
%             fid = fopen(sprintf('texturial_matfiles_vox_warp_gm/%s.aff12.1D',prfix),'w'); % make new local aff
%                 fprintf(fid,'%s\n',newline1);
%                 for i=1:7
%                 fprintf(fid,'%s\n',newline2);
%                 end
%             fclose(fid);
%     
%             %%% %%% %%% %%% GM MASK FIRST
%     
%             infilea = sprintf('texturial_matfiles_vox/%s_gm_Output_Vox_restrict_texmaps_msk.nii.gz',prfix);
%             oufilea = sprintf('texturial_matfiles_vox_warp_gm/%s_gm_mask_deobed.nii.gz',prfix);
%             oufileb = sprintf('texturial_matfiles_vox_warp_gm/%s_gm_mask_reored.nii.gz',prfix);
%             oufilec = sprintf('texturial_matfiles_vox_warp_gm/%s_gm_mask_warped_2mm.nii.gz',prfix);
%     
%             % resliceded it intos the correct orientationalles
%             unix(sprintf('3dWarp -oblique2card -prefix %s -wsinc5 %s', oufilea, infilea)); %-wsinc5
%             unix(sprintf('fslreorient2std %s %s', oufilea, oufileb));        
%             % now we shall do the resamplingue of thine mapses (singlevol)
%             unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 2 -source %s -nwarp "%s/anatQQ_WARP.nii.gz %s/anatQQ.aff12.1D" -prefix %s -ainterp wsinc5',...
%                 warploc,oufileb,warploc,warploc,oufilec));
%     
%             unix(sprintf('rm %s %s',oufilea,oufileb))
%     
%             %%% %%% %%% %%% TEXTURE MAPS NEXT
%     
%             infilea = sprintf('texturial_matfiles_vox/%s_gm_Output_Vox_restrict_texmaps_av.nii.gz',prfix);
%             oufilea = sprintf('texturial_matfiles_vox_warp_gm/%s_tex_maps_deobed.nii.gz',prfix);
%             oufileb = sprintf('texturial_matfiles_vox_warp_gm/%s_tex_maps_reored.nii.gz',prfix);
%             oufilec = sprintf('texturial_matfiles_vox_warp_gm/%s_tex_maps_warped_2mm.nii.gz',prfix);
%     
%             % resliceded it intos the correct orientationalles
%             unix(sprintf('3dWarp -oblique2card -prefix %s -wsinc5 %s', oufilea, infilea)); %-wsinc5
%             unix(sprintf('fslreorient2std %s %s', oufilea, oufileb));        
%             % now we shall do the resamplingue of thine mapses (multivol)
%             unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 2 -source %s -nwarp "%s/anatQQ_WARP.nii.gz texturial_matfiles_vox_warp_gm/%s.aff12.1D" -prefix %s -ainterp wsinc5',...
%                 warploc,oufileb,warploc,prfix,oufilec));
%             
%             unix(sprintf('rm %s %s',oufilea,oufileb))
%         end
%     end
% end
% 
% for i=1:numel(e)
%     i,
%     if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
%             
%         prfix = e(i).name,
%             oufilec = sprintf('texturial_matfiles_vox_warp_gm/%s_gm_mask_warped_2mm.nii.gz',prfix);
%             oufiled = sprintf('texturial_matfiles_vox_warp_gm/%s_gm_mask_warped_2mm_smo6.nii.gz',prfix);
%             unix(sprintf('3dmerge -prefix %s -doall -1blur_fwhm 6 %s',oufiled,oufilec));
% 
%             oufilec = sprintf('texturial_matfiles_vox_warp_gm/%s_tex_maps_warped_2mm.nii.gz',prfix);
%             oufiled = sprintf('texturial_matfiles_vox_warp_gm/%s_tex_maps_warped_2mm_smo6.nii.gz',prfix);
%             unix(sprintf('3dmerge -prefix %s -doall -1blur_fwhm 6 %s',oufiled,oufilec));
%     end
% end


%% SC WARPING
for i=1:numel(e)

    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
            
        prfix = e(i).name,

        if exist(sprintf('texturial_matfiles_vox/%s_sc_Output_Vox_restrict_texmaps_msk.nii.gz',prfix),'file')

            disp('ok!');
        
            warploc = sprintf('Neurocovid/fmri_proc/%s/anat_proc/subpipe_001',prfix);
    
            % create affine for multivolume...
            fid = fopen(sprintf('%s/anatQQ.aff12.1D',warploc),'r'); % read in aff
                newline1 = fgetl(fid);
                newline2 = fgetl(fid);
            fclose(fid);
            fid = fopen(sprintf('texturial_matfiles_vox_warp_sc/%s.aff12.1D',prfix),'w'); % make new local aff
                fprintf(fid,'%s\n',newline1);
                for i=1:7
                fprintf(fid,'%s\n',newline2);
                end
            fclose(fid);

%             %%% %%% %%% %%% GM MASK FIRST (SC1)
%     
%             infilea = sprintf('texturial_matfiles_vox/%s_sc1_Output_Vox_restrict_texmaps_msk.nii.gz',prfix);
%             oufilea = sprintf('texturial_matfiles_vox_warp_sc/%s_sc1_mask_deobed.nii.gz',prfix);
%             oufileb = sprintf('texturial_matfiles_vox_warp_sc/%s_sc1_mask_reored.nii.gz',prfix);
%             oufilec = sprintf('texturial_matfiles_vox_warp_sc/%s_sc1_mask_warped_2mm.nii.gz',prfix);
%     
%             % resliceded it intos the correct orientationalles
%             unix(sprintf('3dWarp -oblique2card -prefix %s -wsinc5 %s', oufilea, infilea)); %-wsinc5
%             unix(sprintf('fslreorient2std %s %s', oufilea, oufileb));        
%             % now we shall do the resamplingue of thine mapses (singlevol)
%             unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 2 -source %s -nwarp "%s/anatQQ_WARP.nii.gz %s/anatQQ.aff12.1D" -prefix %s -ainterp wsinc5',...
%                 warploc,oufileb,warploc,warploc,oufilec));
%     
%             unix(sprintf('rm %s %s',oufilea,oufileb))
%
%             %%% %%% %%% %%% TEXTURE MAPS NEXT (SC1)
%     
%             infilea = sprintf('texturial_matfiles_vox/%s_sc1_Output_Vox_restrict_texmaps_av.nii.gz',prfix);
%             oufilea = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av1_deobed.nii.gz',prfix);
%             oufileb = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av1_reored.nii.gz',prfix);
%             oufilec = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av1_warped_2mm.nii.gz',prfix);
%     
%             % resliceded it intos the correct orientationalles
%             unix(sprintf('3dWarp -oblique2card -prefix %s -wsinc5 %s', oufilea, infilea)); %-wsinc5
%             unix(sprintf('fslreorient2std %s %s', oufilea, oufileb));        
%             % now we shall do the resamplingue of thine mapses (multivol)
%             unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 2 -source %s -nwarp "%s/anatQQ_WARP.nii.gz texturial_matfiles_vox_warp_sc/%s.aff12.1D" -prefix %s -ainterp wsinc5',...
%                 warploc,oufileb,warploc,prfix,oufilec));
%             
%             unix(sprintf('rm %s %s',oufilea,oufileb))
% 
%             %%% %%% %%% %%% TEXTURE MAPS NEXT (SC1-SD)
%     
%             infilea = sprintf('texturial_matfiles_vox/%s_sc1_Output_Vox_restrict_texmaps_sd.nii.gz',prfix);
%             oufilea = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd1_deobed.nii.gz',prfix);
%             oufileb = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd1_reored.nii.gz',prfix);
%             oufilec = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd1_warped_2mm.nii.gz',prfix);
%     
%             % resliceded it intos the correct orientationalles
%             unix(sprintf('3dWarp -oblique2card -prefix %s -wsinc5 %s', oufilea, infilea)); %-wsinc5
%             unix(sprintf('fslreorient2std %s %s', oufilea, oufileb));        
%             % now we shall do the resamplingue of thine mapses (multivol)
%             unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 2 -source %s -nwarp "%s/anatQQ_WARP.nii.gz texturial_matfiles_vox_warp_sc/%s.aff12.1D" -prefix %s -ainterp wsinc5',...
%                 warploc,oufileb,warploc,prfix,oufilec));
%             
%             unix(sprintf('rm %s %s',oufilea,oufileb))


            %%% %%% %%% %%% TEXTURE MAPS NEXT (SC2)
    
            infilea = sprintf('texturial_matfiles_vox/%s_sc2_Output_Vox_restrict_texmaps_av.nii.gz',prfix);
            oufilea = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av2_deobed.nii.gz',prfix);
            oufileb = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av2_reored.nii.gz',prfix);
            oufilec = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av2_warped_2mm.nii.gz',prfix);
    
            % resliceded it intos the correct orientationalles
            unix(sprintf('3dWarp -oblique2card -prefix %s -wsinc5 %s', oufilea, infilea)); %-wsinc5
            unix(sprintf('fslreorient2std %s %s', oufilea, oufileb));        
            % now we shall do the resamplingue of thine mapses (multivol)
            unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 2 -source %s -nwarp "%s/anatQQ_WARP.nii.gz texturial_matfiles_vox_warp_sc/%s.aff12.1D" -prefix %s -ainterp wsinc5',...
                warploc,oufileb,warploc,prfix,oufilec));
            
            unix(sprintf('rm %s %s',oufilea,oufileb))

            %%% %%% %%% %%% TEXTURE MAPS NEXT (SC2-SD)
    
            infilea = sprintf('texturial_matfiles_vox/%s_sc2_Output_Vox_restrict_texmaps_sd.nii.gz',prfix);
            oufilea = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd2_deobed.nii.gz',prfix);
            oufileb = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd2_reored.nii.gz',prfix);
            oufilec = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd2_warped_2mm.nii.gz',prfix);
    
            % resliceded it intos the correct orientationalles
            unix(sprintf('3dWarp -oblique2card -prefix %s -wsinc5 %s', oufilea, infilea)); %-wsinc5
            unix(sprintf('fslreorient2std %s %s', oufilea, oufileb));        
            % now we shall do the resamplingue of thine mapses (multivol)
            unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 2 -source %s -nwarp "%s/anatQQ_WARP.nii.gz texturial_matfiles_vox_warp_sc/%s.aff12.1D" -prefix %s -ainterp wsinc5',...
                warploc,oufileb,warploc,prfix,oufilec));
            
            unix(sprintf('rm %s %s',oufilea,oufileb))


            %%% %%% %%% %%% TEXTURE MAPS NEXT (SC3)
    
            infilea = sprintf('texturial_matfiles_vox/%s_sc3_Output_Vox_restrict_texmaps_av.nii.gz',prfix);
            oufilea = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av3_deobed.nii.gz',prfix);
            oufileb = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av3_reored.nii.gz',prfix);
            oufilec = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av3_warped_2mm.nii.gz',prfix);
    
            % resliceded it intos the correct orientationalles
            unix(sprintf('3dWarp -oblique2card -prefix %s -wsinc5 %s', oufilea, infilea)); %-wsinc5
            unix(sprintf('fslreorient2std %s %s', oufilea, oufileb));        
            % now we shall do the resamplingue of thine mapses (multivol)
            unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 2 -source %s -nwarp "%s/anatQQ_WARP.nii.gz texturial_matfiles_vox_warp_sc/%s.aff12.1D" -prefix %s -ainterp wsinc5',...
                warploc,oufileb,warploc,prfix,oufilec));
            
            unix(sprintf('rm %s %s',oufilea,oufileb))

            %%% %%% %%% %%% TEXTURE MAPS NEXT (SC3-SD)
    
            infilea = sprintf('texturial_matfiles_vox/%s_sc3_Output_Vox_restrict_texmaps_sd.nii.gz',prfix);
            oufilea = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd3_deobed.nii.gz',prfix);
            oufileb = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd3_reored.nii.gz',prfix);
            oufilec = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd3_warped_2mm.nii.gz',prfix);
    
            % resliceded it intos the correct orientationalles
            unix(sprintf('3dWarp -oblique2card -prefix %s -wsinc5 %s', oufilea, infilea)); %-wsinc5
            unix(sprintf('fslreorient2std %s %s', oufilea, oufileb));        
            % now we shall do the resamplingue of thine mapses (multivol)
            unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 2 -source %s -nwarp "%s/anatQQ_WARP.nii.gz texturial_matfiles_vox_warp_sc/%s.aff12.1D" -prefix %s -ainterp wsinc5',...
                warploc,oufileb,warploc,prfix,oufilec));
            
            unix(sprintf('rm %s %s',oufilea,oufileb))

        end
    end
end

for i=1:numel(e)
    i,
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
            
        prfix = e(i).name,

            oufilec = sprintf('texturial_matfiles_vox_warp_sc/%s_sc1_mask_warped_2mm.nii.gz',prfix);
            oufiled = sprintf('texturial_matfiles_vox_warp_sc/%s_sc1_mask_warped_2mm_smo4.nii.gz',prfix);
            unix(sprintf('3dmerge -prefix %s -doall -1blur_fwhm 4 %s',oufiled,oufilec));

            oufilec = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av1_warped_2mm.nii.gz',prfix);
            oufiled = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av1_warped_2mm_smo4.nii.gz',prfix);
            unix(sprintf('3dmerge -prefix %s -doall -1blur_fwhm 4 %s',oufiled,oufilec));

            oufilec = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd1_warped_2mm.nii.gz',prfix);
            oufiled = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd1_warped_2mm_smo4.nii.gz',prfix);
            unix(sprintf('3dmerge -prefix %s -doall -1blur_fwhm 4 %s',oufiled,oufilec));

            oufilec = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av2_warped_2mm.nii.gz',prfix);
            oufiled = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av2_warped_2mm_smo4.nii.gz',prfix);
            unix(sprintf('3dmerge -prefix %s -doall -1blur_fwhm 4 %s',oufiled,oufilec));

            oufilec = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd2_warped_2mm.nii.gz',prfix);
            oufiled = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd2_warped_2mm_smo4.nii.gz',prfix);
            unix(sprintf('3dmerge -prefix %s -doall -1blur_fwhm 4 %s',oufiled,oufilec));

            oufilec = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av3_warped_2mm.nii.gz',prfix);
            oufiled = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_av3_warped_2mm_smo4.nii.gz',prfix);
            unix(sprintf('3dmerge -prefix %s -doall -1blur_fwhm 4 %s',oufiled,oufilec));

            oufilec = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd3_warped_2mm.nii.gz',prfix);
            oufiled = sprintf('texturial_matfiles_vox_warp_sc/%s_tex_maps_sd3_warped_2mm_smo4.nii.gz',prfix);
            unix(sprintf('3dmerge -prefix %s -doall -1blur_fwhm 4 %s',oufiled,oufilec));

    end
end
