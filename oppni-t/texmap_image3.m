function out = texmap_image3( vol, msk, boxwid, dexmax, constr, range_method, step_size, model, b_width, makefigs, subsamp )
% .
% =================================================================
% TEXMAP_IMAGE3: script takes in 3d brain image, sweeps through and gets voxel-by
% voxel measures of local texture values
% ==> for voxelwise mapping only, does a sweep through the image and collects local texture data at back-end 
% =================================================================
%
% out = texmap_image( vol, msk, boxwid, dexmax, constr, range_method, step_size, model, b_width, subsamp )
%
%
% Input:
%
%         vol : 3d image array. Can also trivially handle 2d slices in same format
%         msk : binary ROI mask. must be of same dimensions as "vol". bigger
%               rois may take (slightly) longer to computer features on
%      boxwid : side width of cube centered on each voxel, from which texture is computed
%      dexmax : max. number of dilations to ensure box contains same #voxels as rest
%      constr : constraints on mask growth if numelbnd is specified. Only relevant if numelbnd is non-empty.
%               This is a binary mask, with more coverage than the "msk" volume
%     subsamp : if an integer <N> > 0, applies special "subsampling mode" where texture modelling  
%               is run on only <N> non-overlapping voxel ROIs. used to to
%               fast bandwidth selection by fitting CV curves on a smaller subset of voxels
%     
%      args (range_method, step_size, model, b_width) are passed directly to glca_kde script 
% 
%
% Output:
%
% out.opt_stats : (voxels x 3) array containing median optimal CV error / CV index / bandwidth
% out.CVfull : (bw-values x directions x voxels) array of cross-validation errors
% out.SZfull : (bw-values x directions x voxels) array of bandwidth values
% out.sizecheck : (voxel x 9) array, augmented version of roi_to_glca3 size-checker with preliminary box adjustments included
%        

fprintf('Beginning voxelwise texture analysis ... \n\n');

if nargin<10 || isempty(subsamp)
    subsamp=0;
end

% checks on boxwid variable
if abs(boxwid-floor(boxwid))>eps || boxwid<=1
    error('boxwid must be odd integer, greater than 1');
elseif mod(boxwid,2)==0
    error('boxwid must be odd integer');
elseif size(vol,3)>1 && boxwid > 1.0*min(size(vol))
    error('boxwid is *very* big relative to 3d image size (spans >100% of smallest dim)');
elseif size(vol,3)==1 && boxwid > 1.0*min([size(vol,1) size(vol,2)])
    error('boxwid is *very* big relative to 2d image size (spans >100% of smallest dim)');
end
% boxwidth, converted into number of steps from center / bidirectional
d0 = (boxwid-1)/2;

% now get linear indices of voxels in mask
vix = find( msk>0 );
[ia,ja,ka]=ind2sub( size(msk), vix );

if subsamp==0 % full sweep - texture estimation
    
    notepts = unique(round(linspace(1,numel(vix),100)));
    
    tic;
    for v=1:numel(vix) 
    
        if sum( v==notepts )>0
           toc,
           fprintf('   completed %.01f pct of voxels (%u/%u)...\n', round(100*v/numel(vix)), v, numel(vix) );
        end
        
        db = 0;
        sz_flag = 0;
        while sz_flag==0 && db < dexmax
    
                ijk1 = [ia(v) ja(v) ka(v)]-(d0+db); % start idx
                ijk1 = max( [ijk1; [1 1 1]],[],1); % push back to 1
                ijk2 = [ia(v) ja(v) ka(v)]+(d0+db); % end idx
                ijk2 = min( [ijk2; [size(vol)]],[],1); % push back to maxdim
                % constraint-box
                cstbox = constr( ijk1(1):ijk2(1), ijk1(2):ijk2(2), ijk1(3):ijk2(3) );
                if db==0
                    init_numvox = sum(cstbox(:));
                end
    
            if sum(cstbox(:)) >= boxwid^3
                sz_flag = 1;
            else
                db = db+1;
            end
        end
        finl_numvox = sum(cstbox(:));
    
        % take box chunk - should also handle 2d ok
        volbox = vol( ijk1(1):ijk2(1), ijk1(2):ijk2(2), ijk1(3):ijk2(3) );
        mskbox = ones(size(volbox));
        
        % run texture analysis
        outtmp = roi_to_glca3( volbox, mskbox, boxwid^3, cstbox, range_method, step_size, model, b_width, makefigs );
    
        out.metrics_av(v,:) = outtmp.metrics_av;
        out.metrics_sd(v,:) = outtmp.metrics_sd;
        out.bscalset(1)     = outtmp.bscalset(1);
        out.Bvl_for_fitt    = outtmp.Bvl_for_fitt;
        % augmenting size checker with preliminary box adjustments
        % [upper & low bounds][init/final/pct-cht --first resize][init/final/pct-cht --last resize][warn-]
        out.sizecheck(v,:)  = [outtmp.sizecheck(1:2),[init_numvox, finl_numvox, 100*(finl_numvox-init_numvox)/init_numvox],outtmp.sizecheck(3:6)];
    end

else

    if subsamp>0.25*numel(vix)
        error('subsampling number is big compared to overall sample size');
    end
    vix_list = randperm(numel(vix));

    tic;
    kq=0;
    ijka_coll = [99999 99999 99999]; % collected coordinates
    for v=1:subsamp
    
        fprintf('   completed %.01f pct of voxels (%u/%u)...\n', round(100*v/subsamp), v, subsamp );

        kq=kq+1; % step to next vox in list

        % ====loop until you find a non-overlapping centroid
        ijka_cand = [ia(vix_list(kq)) ja(vix_list(kq)) ka(vix_list(kq))]; %-candidate
        while max(sum( abs(bsxfun(@minus,ijka_coll,ijka_cand)) <= boxwid,2 ))==3
            kq=kq+1;
            kq,
            if kq>numel(vix_list)
               error('cannot find enough spaced rois to sample -- collect less!')
            end
            ijka_cand = [ia(vix_list(kq)) ja(vix_list(kq)) ka(vix_list(kq))]; %-candidate retry on next
        end
        % ====loop finished
        
        % ====loop until you find a box that can reach the correct size
        db = 0; %-candidate
        sz_flag = 0;
        while sz_flag==0 && db < dexmax
    
                ijk1 = [ia(vix_list(kq)) ja(vix_list(kq)) ka(vix_list(kq))]-(d0+db); % start idx
                ijk1 = max( [ijk1; [1 1 1]],[],1); % push back to 1
                ijk2 = [ia(vix_list(kq)) ja(vix_list(kq)) ka(vix_list(kq))]+(d0+db); % end idx
                ijk2 = min( [ijk2; [size(vol)]],[],1); % push back to maxdim
                % constraint-box
                cstbox = constr( ijk1(1):ijk2(1), ijk1(2):ijk2(2), ijk1(3):ijk2(3) );
                if db==0
                    init_numvox = sum(cstbox(:));
                end
    
            if sum(cstbox(:)) >= boxwid^3
                sz_flag = 1;
                finl_numvox = sum(cstbox(:));
            else
                db = db+1;
            end
        end

        ijka_coll = [ijka_coll; [ia(vix_list(kq)) ja(vix_list(kq)) ka(vix_list(kq))]];

        % take box chunk - should also handle 2d ok
        volbox = vol( ijk1(1):ijk2(1), ijk1(2):ijk2(2), ijk1(3):ijk2(3) );
        mskbox = ones(size(volbox));

        % run texture analysis
        outtmp = roi_to_glca3( volbox, mskbox, boxwid^3, cstbox, range_method, step_size, model, b_width, makefigs );
    
        out.opt_stats(v,:) = [median(outtmp.Copt_err) median(outtmp.Copt_idx) median(outtmp.Copt_bvl)];
        out.CVfull(:,:,v) = outtmp.CVfull;
        out.SZfull(:,:,v) = outtmp.SZfull;
        % augmenting size checker with preliminary box adjustments
        out.sizecheck(v,:)  = [outtmp.sizecheck(1:2),[init_numvox, finl_numvox, 100*(finl_numvox-init_numvox)/init_numvox],outtmp.sizecheck(3:6)];
    end
    
end

fprintf('\nFinished!\n');


% % % . temp plotting stuff
% % for z = 120;
% % vol_slc = vol(:,:,z);
% % msk_slc = msk(:,:,z);
% % edg = edge(msk_slc);
% % bdg = zeros(size(edg));
% % bdg((100-3-1):(100+3+1),112-3-1) = 1;
% % bdg((100-3-1):(100+3+1),112+3+1) = 1;
% % bdg(100-3-1,(112-3-1):(112+3+1)) = 1;
% % bdg(100+3+1,(112-3-1):(112+3+1)) = 1;
% % rron = cat(3,ones(size(vol_slc)),zeros(size(vol_slc)),zeros(size(vol_slc)));
% % bron = cat(3,zeros(size(vol_slc)),zeros(size(vol_slc)),ones(size(vol_slc)));
% % %figure,imagesc( vol_slc ); colormap gray; hold on; h=imagesc( gron); set(h,'AlphaData',0.2*msk_slc);
% % figure,imagesc( vol_slc ); colormap gray; hold on; 
% % h=imagesc( rron); set(h,'AlphaData',1*edg);
% % l=imagesc( bron); set(l,'AlphaData',1*bdg);
% % end
% % ixvv = find( ia==100 & ja==112 & ka==120)
% % %
% % figure,imagesc( volbox(:,:,4)); colormap gray;

