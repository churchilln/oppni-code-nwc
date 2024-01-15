function out = spike_estimator( volname, maskname, mpename, outprefix, WSZ, BIDIR )
%
% just estimates the spikes, doesn't remove them yet
%
%

% load and mean-center the mpes
mpemat = load(mpename);
mpemat = bsxfun(@minus,mpemat,mean(mpemat));

% load and mean-center the fmri data
VV = load_untouch_niiz(volname); % load 4D fMRI dataset
MM = load_untouch_niiz(maskname);% load the binary brain mask
% convert mask into 3D matrix volume
mask = double(MM.img);
% convert 4D fMRI data into (voxels x time) matrix
epimat = nifti_to_mat( VV,MM ); 
epimat = bsxfun(@minus,epimat,mean(epimat,2));

% matrix dimensions
[Ntime1 Nvar] = size(mpemat);
[Nvox Ntime2] = size(epimat);
if( Ntime1 == Ntime2 ) Ntime=Ntime1; clear Ntime1 Ntime2;
else error('different timepoints in fmri and mpe');
end

%% for mpes...

[u,l,v]=svd(mpemat','econ');
QQ = (v(:,1:end-1)*l(1:end-1,1:end-1))';

dist_mot = zeros(Ntime,1);
for(t=1:Ntime)
    if(BIDIR>0) % ==> bidirectional window
        window = (t-WSZ):(t+WSZ);
    else % ==> monodirectional window
        window = (t-WSZ):t;
    end
    % trim overrun, and self from window
    window(window<1)     = [];
    window(window>Ntime) = [];
    window(window==t)    = [];
    % get distance score from medioid
    if(~isempty(window))
    dist_mot(t) = sum( (QQ(:,t) - median(QQ(:,window),2)).^2 );
    end
end
% convert to prob.
xtmp   = dist_mot+eps; xtmp=xtmp./max(xtmp);
par_ab = gamfit( xtmp );
prob_mot = 1-gamcdf( xtmp, par_ab(1), par_ab(2) );

%% for bold(whole-vol)...

[u,l,v]=svd(epimat,'econ');
QQ = (v(:,1:end-1)*l(1:end-1,1:end-1))';

dist_vol = zeros(Ntime,1);
for(t=1:Ntime)
    if(BIDIR>0) % ==> bidirectional window
        window = (t-WSZ):(t+WSZ);
    else % ==> monodirectional window
        window = (t-WSZ):t;
    end
    % trim overrun, and self from window
    window(window<1)     = [];
    window(window>Ntime) = [];
    window(window==t)    = [];
    % get distance score from medioid
    if(~isempty(window))
    dist_vol(t) = sum( (QQ(:,t) - median(QQ(:,window),2)).^2 );
    end
end
% convert to prob.
xtmp   = dist_vol+eps; xtmp=xtmp./max(xtmp);
par_ab = gamfit( xtmp );
prob_vol = 1-gamcdf( xtmp, par_ab(1), par_ab(2) );

%% for bold(slice-based)...

% get 3D volume dimensions
[Nx Ny Nz] = size( mask );
% get index-vector, where each voxel is labelled with its slice# (z=1...Nz)
ixvol = ones( [Nx Ny Nz] );
for(z=1:Nz) ixvol(:,:,z) = ixvol(:,:,z) .* z; end
ixvect = ixvol(mask>0);

dist_slc = zeros(Ntime,Nz);
for(z=1:Nz)
    
    % check for signal in slice
    fullsum = sum(sum(abs(epimat(ixvect==z,:))));
    
    % requires that there be signal in this slice
    if( (fullsum~=0) && isfinite( fullsum ) )
        [u,l,v]=svd(epimat(ixvect==z,:),'econ');
        QQ = (v(:,1:end-1)*l(1:end-1,1:end-1))';
        
        for(t=1:Ntime)
            if(BIDIR>0) % ==> bidirectional window
                window = (t-WSZ):(t+WSZ);
            else % ==> monodirectional window
                window = (t-WSZ):t;
            end
            % trim overrun, and self from window
            window(window<1)     = [];
            window(window>Ntime) = [];
            window(window==t)    = [];
            % get distance score from medioid
            if(~isempty(window))
            dist_slc(t,z) = sum( (QQ(:,t) - median(QQ(:,window),2)).^2 );
            end
        end
        % convert to prob.
        xtmp   = dist_slc(:,z)+eps; xtmp=xtmp./max(xtmp);
        par_ab = gamfit( xtmp );
        prob_slc(:,z) = 1-gamcdf( xtmp, par_ab(1), par_ab(2) );

    else
        dist_slc(:,z) = NaN*ones(Ntime,1);
        prob_slc(:,z) =     ones(Ntime,1);
    end
end

% time x wind x type
out.outl_mot = double( prob_mot < 0.05 );
out.outl_vol = double( prob_vol < 0.05 );
out.outl_slc = double( prob_slc < 0.05 );
out.outl_vol_mot = double( out.outl_vol & ( [0; out.outl_mot(1:end-1)] | [out.outl_mot(2:end); 0] ) );
for(z=1:Nz)
out.outl_slc_mot(:,z) = double( out.outl_slc(:,z) & ( [0; out.outl_mot(1:end-1)] | [out.outl_mot(2:end); 0] ) );
end
% also store actual disp traces
out.ssd_mot = dist_mot;
out.ssd_vol = dist_vol;
out.ssd_slc = dist_slc;

% figure,
% subplot(3,2,1); plot( dist_mot );
% subplot(3,2,3); plot( dist_vol );
% subplot(3,2,5); plot( dist_slc );
% subplot(2,2,2); imagesc( dist_slc );
% 
% figure,
% subplot(3,2,1); bar( out.outl_mot );
% subplot(3,2,3); bar( out.outl_vol );
% subplot(3,2,5); bar( out.outl_slc );
% subplot(2,2,2); imagesc( out.outl_slc );
% 
% figure,
% subplot(3,2,1); bar( out.outl_mot );
% subplot(3,2,3); bar( out.outl_vol_mot );
% subplot(3,2,5); bar( out.outl_slc_mot );
% subplot(2,2,2); imagesc( out.outl_slc_mot );
