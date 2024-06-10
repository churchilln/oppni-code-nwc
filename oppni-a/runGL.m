function runGL(bg, overlay, imname, config)
%
%  runGL(bg, overlay, imname, config)
%
%  
%  imname = {prefix of figure}
%  config = 'anat_mask_unwarp', 'func_mask_unwarp', 'anat_warp','func_warp'
%
%


unix(sprintf('rm -r QCFigs_%s',config));
unix(sprintf('mkdir QCFigs_%s',config));

%launch MRIcroGL and display an image
% bg: name of background image to load
% overlay: name of overlay image
glExe = '/Applications/MRIcroGL.app/Contents/MacOS/MRIcroGL';
if ~exist(glExe,'file'), error('Unable to find %s', glExe); end
if ~exist('bg','var'), bg = 'spm152'; end
if ~exist('overlay','var'), overlay = 'spmMotor'; end

%create script
fnm = fullfile(fileparts(mfilename),'script.py');
fileID = fopen(fnm,'w');

fprintf(fileID, 'import gl\n');
fprintf(fileID, 'gl.resetdefaults()\n');

for i=1:numel(bg)

    pngout = sprintf('QCFigs_%s/%s.png',config,imname{i});

    fprintf(fileID, sprintf('gl.loadimage("%s")\n', bg{i}));
    fprintf(fileID, sprintf('gl.overlayload("%s")\n', overlay{i}));

    %---> display

    if strcmpi(config,'anat_mask_unwarp')
        
        if contains(overlay{i},'.nii.gz')
            ovlfile = '__opptmp_rungl_anat.nii';
            unix(sprintf('cp %s %s.gz',overlay{i},ovlfile));
            unix(sprintf('gunzip %s.gz',ovlfile));
        elseif contains(overlay{i},'.nii')
            ovlfile = overlay{i};
        else
            error('unrecognized extension!')
        end

        V=load_untouch_niiz(ovlfile);
        lo = V.hdr.hist.qoffset_z;
        sz = V.hdr.dime.pixdim(4);
        nm = V.hdr.dime.dim(4);
        slcsum = squeeze(sum(sum(abs(double(V.img)),1),2));
        ixbnd = find( slcsum > 0.01*max(slcsum) );
        ixbnd = ixbnd([1 end]);
        ixbnd = (ixbnd*sz) + lo;
        slcvals = round( linspace( ixbnd(1), ixbnd(2), 21 ) );

        if(contains(overlay{i},'.gz'))
            unix(sprintf('rm %s',ovlfile));
        end
        
        fprintf(fileID, sprintf('gl.opacity(1,50)\n'));
        fprintf(fileID, sprintf('gl.mosaic("A H 0.2 V 0.2 %d %d %d %d %d %d %d; %d %d %d %d %d %d %d; %d %d %d %d %d %d %d")\n',slcvals));
        fprintf(fileID, sprintf('gl.savebmp(''%s'')\n',pngout));

    elseif strcmpi(config,'func_mask_unwarp')
        overlay{i},
        if contains(overlay{i},'.nii.gz')
            ovlfile = '__opptmp_rungl_anat.nii';
            unix(sprintf('cp %s %s.gz',overlay{i},ovlfile));
            unix(sprintf('gunzip %s.gz',ovlfile));
        elseif contains(overlay{i},'.nii')
            ovlfile = overlay{i};
        else
            error('unrecognized extension!')
        end

        V=load_untouch_niiz(ovlfile);
        lo = V.hdr.hist.qoffset_z;
        sz = V.hdr.dime.pixdim(4);
        nm = V.hdr.dime.dim(4);
        slcsum = squeeze(sum(sum(abs(double(V.img)),1),2));
        ixbnd = find( slcsum > 0.01*max(slcsum) );
        ixbnd = ixbnd([1 end]);
        ixbnd = (ixbnd*sz) + lo;
        slcvals = round( linspace( ixbnd(1), ixbnd(2), 21 ) );

        if(contains(overlay{i},'.gz'))
            unix(sprintf('rm %s',ovlfile));
        end
        
        fprintf(fileID, sprintf('gl.opacity(1,50)\n'));
        fprintf(fileID, sprintf('gl.mosaic("A H 0.2 V 0.2 %d %d %d %d %d %d %d; %d %d %d %d %d %d %d; %d %d %d %d %d %d %d")\n',slcvals));
        fprintf(fileID, sprintf('gl.savebmp(''%s'')\n',pngout));

    elseif strcmpi(config,'anat_warp')

        if(contains(bg{i},'.nii.gz'))
            ovlfile = '__opptmp_rungl_anat.nii';
            unix(sprintf('cp %s %s.gz',bg{i},ovlfile));
            unix(sprintf('gunzip %s.gz',ovlfile));
        elseif contains(bg{i},'.nii')
            ovlfile = bg{i};
        else
            error('unrecognized extension!')
        end

        V=load_untouch_niiz(ovlfile);
        bnd = double( prctile(V.img(V.img>eps),[2.5 97.5]));

        if(contains(bg{i},'.gz'))
            unix(sprintf('rm %s',ovlfile));
        end
        
        fprintf(fileID, sprintf('gl.minmax(0,%d,%d)\n',round(bnd)));
        fprintf(fileID, sprintf('gl.opacity(1,33)\n'));
        fprintf(fileID, sprintf('gl.mosaic("S -50 C -65 A -6; S -8 C -20 A 13; S 30 C 54 A 58")\n'));
        fprintf(fileID, sprintf('gl.savebmp(''%s'')\n',pngout));

    elseif strcmpi(config,'func_warp')

        if(contains(bg{i},'.gz'))
            ovlfile = '__opptmp_rungl_anat.nii';
            unix(sprintf('cp %s %s.gz',bg{i},ovlfile));
            unix(sprintf('gunzip %s.gz',ovlfile));
        elseif contains(bg{i},'.nii')
            ovlfile = bg{i};
        else
            error('unrecognized extension!')
        end

        V=load_untouch_niiz(ovlfile);
        bnd = double( prctile(V.img(V.img>eps),[2.5 97.5]));

        if(contains(bg{i},'.gz'))
            unix(sprintf('rm %s',ovlfile));
        end
        
        fprintf(fileID, sprintf('gl.minmax(0,%d,%d)\n',round(bnd)));
        fprintf(fileID, sprintf('gl.opacity(1,33)\n'));
%         fprintf(fileID, sprintf('gl.mosaic("A -30 -19 -8; 3 14 25; 35 46 57; 68 79 90")\n'));
        fprintf(fileID, sprintf('gl.mosaic("S -50 C -65 A -6; S -8 C -20 A 13; S 30 C 54 A 58")\n'));
        fprintf(fileID, sprintf('gl.savebmp(''%s'')\n',pngout));
    end
end

%fprintf(fileID, sprintf('gl.quit()'));
fclose(fileID);

%run script
cmd = [glExe,' "', fnm ,'" &']
system(cmd);

return;

% % % catpath = 'AFNI.MNI152_2009_template_SSW.3mm';
% % % outpath = '/Users/tomschweizer/Documents/MATLAB/Neurocovid/fmri_proc';
% % % 
% % % Va = load_untouch_niiz([outpath,'/_group_level/brain_maps/',catpath,'/anat_brain_grp.nii']);
% % % Vb = load_untouch_niiz([outpath,'/_group_level/masks/',catpath,'/anat_brain_mask_grp.nii']);
% % % V  = Vb;
% % % vsm = double(Va.img);
% % % vsi = 1 - (vsm - min(vsm(:)))./(max(vsm(:))-min(vsm(:)));
% % % vsi  = vsi .* double( Vb.img );
% % % V.img = double( vsi > 0.4);% prctile(vsi(vsi>0),75));
% % % save_untouch_niiz(V,[outpath,'/_group_level/brain_maps/',catpath,'/anat_sulcal_grp.nii'])

% % catpath = 'AFNI.MNI152_2009_template_SSW.3mm';
% % outpath = '/Users/tomschweizer/Documents/MATLAB/Neurocovid/fmri_proc';
% % 
% % Va = load_untouch_niiz([outpath,'/_group_level/brain_maps/',catpath,'/func_tAV_grp.nii']);
% % Vb = load_untouch_niiz([outpath,'/_group_level/masks/',catpath,'/func_brain_mask_grp.nii']);
% % V  = Vb;
% % vsm = double(Va.img);
% % vsi = (vsm - min(vsm(:)))./(max(vsm(:))-min(vsm(:)));
% % vsi  = vsi .* double( Vb.img );
% % V.img = double( vsi > 0.55);% prctile(vsi(vsi>0),75));
% % save_untouch_niiz(V,[outpath,'/_group_level/brain_maps/',catpath,'/func_sulcal_grp.nii'])
