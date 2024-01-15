function texmap_wrapper3( image, mask, masklab, isvoxmap, outfile, BOXWID, NBND, range_method, step_size, model, b_width, makefigs, mergmode )
% .
% =================================================================
% TEXMAP_WRAPPER: wrapper script that takes in nifti volume, mask, and some other information
% and generates texture maps : as both 4d nifti volume set and vectorized matlab objects
% =================================================================
%
% Syntax:
%          texmap_wrapper3( image, mask, masklab, isvoxmap, outfile, BOXWID, NBND, range_method, step_size, model, b_width, makefigs, mergmode )
%
% Input:
%          Image : string specifying input image file
%          Mask : string specifying input brain mask (or parcellation)
%          Masklab : if Mask is a parcellation, you can input integer list of parcels to include as mask
%          Isvoxmap : binary, indicating whether to do voxel-wise mapping
%          Outfile : string specifying the path/name of outputs
%          BOXWID : if voxel-wise mapping, width of "box" for texture estimation (must be an odd integer)
%          NBND : bounds on desired roi mask size.
%          mergmode : if roi-based mapping, merge masklab parcels into single mask instead of handling separately 
%   
%          (args range_method, step_size, model, b_width, makefigs) are passed to glca_kde
%
%          Note:
%                > if masklab =[], isvoxmap=1 --> voxelwise over all mask areas>0                           (code: Vox_allmask)
%                > if masklab =[], isvoxmap=0 --> global over all mask>0                                    (code: Roi_allmask)
%                > if masklab~=[], isvoxmap=1 --> voxelwise, restricted to union of listed parc labels      (code: Vox_restrict) 
%                > if masklab~=[], isvoxmap=0 --> global, computed separately within each listed parc label (code: Roi_byparc) 
%
% Output:
%          a matfile containing texture data and, in case of voxelwise
%          mapping, nifti files with maps
%
%
%          sizecheck:
%                       [(lower bound) (upper bound), (initial size) (final size - coarse adjust) (percent change - coarse adjust), (initial size, post-coarse) (final size - fine adjust) (percent change - fine adjust), (warning)]  
%
%


if nargin < 12
    makefigs = 0;
end
if nargin < 13
    mergmode = 0;
end

if isvoxmap==0
    if isempty(masklab)
        eu='Roi_allmask';
    else
        eu='Roi_byparc';
    end
elseif isvoxmap==1
    if isempty(masklab)
        eu='Vox_allmask';
    else
        eu='Vox_restrict';
    end
else
    error('unrecognized voxel-map option')
end

filostr = sprintf('%s_Output_%s',outfile,eu);

if ~exist( sprintf('%s.mat',filostr), 'file')

    V = load_untouch_niiz(image);
    A = load_untouch_niiz(mask);
    
    volimag = double(V.img);
    gm_mask = zeros(size(A.img)); % initialize the mask

    if isvoxmap==0 % roi based mapping

        if isempty(masklab) || mergmode>0
            if isnumeric(b_width) && numel(b_width)~=1
                error('global roi mapping / roi mapping in mergmode can only take single bw estimate');
            end
            if isempty(masklab)
                gm_mask = double(A.img>eps);
            elseif mergmode>0
                gm_mask = 0;
                for v=1:numel(masklab)
                    gm_mask = gm_mask + double( A.img==masklab(v) );
                end
                gm_mask = double(gm_mask>0);
            end

            v=1;
                'running',
                [v numel(masklab)],
                %>> gm_mask was pre-specified
                sum(gm_mask(:)),
                % run texture analysis
                if ~isnumeric(b_width)
                outtmp =roi_to_glca3( volimag, gm_mask, NBND, gm_mask, range_method, step_size, model, b_width, makefigs );        
                else
                outtmp =roi_to_glca3( volimag, gm_mask, NBND, gm_mask, range_method, step_size, model, b_width(v), makefigs );        
                end
                % put into array for output
                if ~isnumeric(b_width) && strcmpi(b_width,'estim')
                    disp(' ...no metrics to save...')
                else
                    out.metrics_av(v,:) = outtmp.metrics_av;
                    out.metrics_sd(v,:) = outtmp.metrics_sd;
                    out.bscalset(v,:)   = outtmp.bscalset;
                    if strcmpi(model,'KDE')
                        out.TFopt(v,1) = outtmp.TFopt;
                    end
                end
                if ~isnumeric(b_width) && (strcmpi(b_width,'optim') || strcmpi(b_width,'estim'))
                    out.opt_stats(v,:)  = [median(outtmp.Copt_err) median(outtmp.Copt_idx) median(outtmp.Copt_bvl)];
                    out.CVfull(:,:,v) = outtmp.CVfull;
                    out.SZfull(:,:,v) = outtmp.SZfull;
                    if strcmpi(model,'KDE')
                        out.TFset(:,:,v) = outtmp.TFset;
                    end
                else
                    disp(' ...no optimization stats tot save...')
                end
                out.Bvl_for_fitt(v,1) = outtmp.Bvl_for_fitt;
                % augmented to match voxel approach
                out.sizecheck(v,:)  = [outtmp.sizecheck(1:2), [NaN NaN NaN], outtmp.sizecheck(3:6)];
            v=1;
        else
            if isnumeric(b_width) && numel(b_width)~=1 && numel(b_width)~=numel(masklab)
                error('regional roi mapping can only take single bw estimate / matching #labels');
            elseif isnumeric(b_width) && numel(b_width)==1 % if single value, copy to match roi label count
                b_width = repmat( b_width, numel(masklab), 1);
            end
            %gm_sum=0;
            if size(NBND,1)~=numel(masklab)
                error('boundary array doesnt match number of mask labels');
            end

            for v=1:numel(masklab)
                'running',
                [v numel(masklab)],
                gm_mask = double( A.img==masklab(v) );
                sum(gm_mask(:)),
                % run texture analysis
                if ~isnumeric(b_width)
                outtmp =roi_to_glca3( volimag, gm_mask, NBND(v,:), gm_mask, range_method, step_size, model, b_width, makefigs );        
                else
                outtmp =roi_to_glca3( volimag, gm_mask, NBND(v,:), gm_mask, range_method, step_size, model, b_width(v), makefigs );        
                end
                % put into array for output
                if ~isnumeric(b_width) && strcmpi(b_width,'estim')
                    disp(' ...no metrics to save...')
                else
                    out.metrics_av(v,:) = outtmp.metrics_av;
                    out.metrics_sd(v,:) = outtmp.metrics_sd;
                    out.bscalset(v,:)   = outtmp.bscalset;
                    if strcmpi(model,'KDE')
                        out.TFopt(v,1) = outtmp.TFopt;
                    end
                end
                if ~isnumeric(b_width) && (strcmpi(b_width,'optim') || strcmpi(b_width,'estim'))
                    out.opt_stats(v,:)  = [median(outtmp.Copt_err) median(outtmp.Copt_idx) median(outtmp.Copt_bvl)];
                    out.CVfull(:,:,v) = outtmp.CVfull;
                    out.SZfull(:,:,v) = outtmp.SZfull;
                    if strcmpi(model,'KDE')
                        out.TFset(:,:,v) = outtmp.TFset;
                    end
                else
                    disp(' ...no optimization stats tot save...')
                end
                out.Bvl_for_fitt(v,1) = outtmp.Bvl_for_fitt;
                % augmented to match voxel approach
                out.sizecheck(v,:)  = [outtmp.sizecheck(1:2), [NaN NaN NaN], outtmp.sizecheck(3:6)];
            end
        end

    elseif isvoxmap==1 % voxelwise mapping

        % --currently only supports HDE-based estimation for voxelwise
        if ~strcmpi(model,'HDE')
            error('only the HDE-based estimator model currently supported for voxelwise analysis. Kernel methods are too slow!');
        end

        % --currently only supports scalar BW ... spatial map maybe in future
        if isnumeric(b_width) && numel(b_width)~=1
            error('voxelwise texture analysis only supports global BW estimator for now')
        end

        if isempty(masklab)
            gm_mask = double(A.img>eps);
        else
            for i=1:numel(masklab)
                gm_mask = gm_mask + double( A.img==masklab(i) );
            end
            gm_mask(gm_mask>1)=1;
        end
        % this is in predef mode:
        if ~isnumeric(b_width)
            if strcmpi(b_width,'estim')
                % exports the subsetted voxelwise list for estim-mode: opt_stats, CVfull, SZfull  
                out = texmap_image3( volimag, gm_mask, BOXWID, BOXWID, gm_mask, range_method, step_size, model, b_width, makefigs, 20 );
            elseif strcmpi(b_width,'optim')
                error('full optimization not supported for voxelwise mapping!')
            end
        else
            % exports the voxelwise list for prespecced BW:  metrics_av/metrics_sd/bscalset/sizecheck
            
            % [upper & low bounds][init/final/pct-cht --first resize][init/final/pct-cht --last resize][warn-]
            out = texmap_image3( volimag, gm_mask, BOXWID, BOXWID, gm_mask, range_method, step_size, model, b_width, makefigs, 0 );
        end
    else
        error('unrecognized "isvoxmap" setting');
    end
    
    save(sprintf('%s.mat',filostr),'out');
    unix(sprintf('rm __texmap_intermed/image_%s.nii.gz',filostr));
    unix(sprintf('rm __texmap_intermed/mask_%s.nii.gz',filostr));

else
    fprintf('matfile already exists for "%s". Delete if you want to rerun.\n',filostr);
end

%% pushing to niftis -- voxelwise mapping only!

if isvoxmap==1 && isnumeric(b_width) && (~exist( sprintf('%s_texmaps_av.nii.gz',filostr), 'file') || ~exist( sprintf('%s_texmaps_sd.nii.gz',filostr), 'file'))
    disp('running nii vol reconstruction...')

    V = load_untouch_niiz(image);
    A = load_untouch_niiz(mask);
    
    volimag = double(V.img);
    gm_mask = zeros(size(A.img)); % initialize the mask

    % modify mask per label set
    if isempty(masklab)
        gm_mask = double(A.img>eps);
    else
        for i=1:numel(masklab)
            gm_mask = gm_mask + double( A.img==masklab(i) );
        end
        gm_mask(gm_mask>1)=1;
    end

    load(sprintf('%s.mat',filostr));
    param_list = {'ENE','ENT','CON','HOM','COR','SHA','PRO'};

    % making 7 "AVG" map-files
    TMPVOL=zeros( [size(gm_mask) 7]);
    for ip=1:7
        tmp = gm_mask;
        tmp(tmp>0) = out.metrics_av(:,ip);
        TMPVOL(:,:,:,ip) = tmp;
    end
    % port to nifti
    nii=V;
    nii.img = TMPVOL;
    nii.hdr.dime.datatype = 16;
    nii.hdr.hist = V.hdr.hist;
    nii.hdr.dime.dim(1) = 4;
    nii.hdr.dime.dim(5) = size(TMPVOL,4);
    save_untouch_niiz(nii,sprintf('%s_texmaps_av.nii.gz',filostr)); 

    % making 7 "STD" map-files
    TMPVOL=zeros( [size(gm_mask) 7]);
    for ip=1:7
        tmp = gm_mask;
        tmp(tmp>0) = out.metrics_sd(:,ip);
        TMPVOL(:,:,:,ip) = tmp;
    end
    % port to nifti
    nii=V;
    nii.img = TMPVOL;
    nii.hdr.dime.datatype = 16;
    nii.hdr.hist = V.hdr.hist;
    nii.hdr.dime.dim(1) = 4;
    nii.hdr.dime.dim(5) = size(TMPVOL,4);
    save_untouch_niiz(nii,sprintf('%s_texmaps_sd.nii.gz',filostr)); 

    % and corresponding mask reference
    TMPVOL = gm_mask;
    nii=V;
    nii.img = TMPVOL;
    nii.hdr.dime.datatype = 16;
    nii.hdr.hist = V.hdr.hist;
    nii.hdr.dime.dim(1) = 4;
    nii.hdr.dime.dim(5) = size(TMPVOL,4);
    save_untouch_niiz(nii,sprintf('%s_texmaps_msk.nii.gz',filostr)); 
else
    fprintf('nifti files already exist for "%s". Delete if you want to rerun.\n',filostr);
end
