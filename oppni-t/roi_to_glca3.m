function out = roi_to_glca3( vol, msk, numelbnd, constr, range_method, step_size, model, b_width, makefigs )
% .
% =================================================================
% ROI_GLCA3: script takes in 3d image array, binary ROI mask (same size),
% and some constraints on ROI growing/shrinking. Collects list of texture measures within 
% size-adjusted/masked-out ROI
% =================================================================
%
% Syntax:
%           out = roi_to_glca( vol, msk, numelbnd, constr, range_method, step_size, model, b_width, makefigs )
%
% input:
%
%         vol : 3d image array. Can also trivially handle 2d slices in same format
%         msk : binary ROI mask. must be of same dimensions as "vol". bigger
%               rois may take (slightly) longer to computer features on
%    numelbnd : desired roi mask size. If ~isempty(numelbnd), and your input mask is too big/too small,
%               it will be first dilated to meet/exceed numelbnd (if needed),
%               then will prune voxels, starting at those most distant from the center of mass, until 
%               exactly matches numelbnd (if needed). If numelbnd = [], leaves input mask unchanged
%      constr : constraints on mask growth if numelbnd is specified. Only relevant if numelbnd is non-empty.
%               This is a binary mask, with more coverage than the "msk" volume
%
%        args (range_method, step_size, model, b_width, makefigs) are passed directly to glca_kde script 
%
% output:
%	   out.metrics_av : average of texture metrics, glca script outputs
%	   out.metrics_sd : stdev of texture metrics, glca script outputs
%       out.sizecheck : vector of information related to ROI resizing.
%                       includes:
%                       [[numelbnd] init_numvox, finl_numvox, pct_numvox, dwarn]
%
%                       where init_numvox is the number of input voxels,
%                       finl_numvox is the number of voxels sent to
%                       analysis after dilating/trimming, pct_numvox is the
%                       percent change from init to finl, and dwarn is a
%                       binary flag -- warns if algorithm fails to make mask large enough to exceed numelbnd
%                       within pre-specified number of dilations (currently set at 20)

 maxdil = 20; % hard-coded limit --

% possible choices
% -unmodified vol + mask
% -vol+mask with size bound
% -vol+mask with size bound & constraints on growth

if nargin<3 % --== no changing of roi
    numelbnd=[];
end
if nargin<4 || isempty(constr) % --== no constraints on resizing
   constr = 1; % becomes simple scalar instead
end
if nargin<9 || isempty(makefigs)
    makefigs=0;
end
        
% dimensions of mask array
[nx,ny,nz]=size(msk);
% bounded by constrained area
msk = msk .* constr;
% record for later check
init_numvox = sum(msk(:));

%% (1) resize ROI for fixed volume estimation > texture on consistently sized rois 

% --
% tic;
if ~isempty(numelbnd) % if we are bounding roi size...
    %disp('resizing');
    % record roi center of mass after bounding
    d1=sum(sum(msk,2),3); mcom(1) = sum(d1(:).*(1:size(msk,1))')./sum(d1(:));
    d1=sum(sum(msk,1),3); mcom(2) = sum(d1(:).*(1:size(msk,2))')./sum(d1(:));
    d1=sum(sum(msk,1),2); mcom(3) = sum(d1(:).*(1:size(msk,3))')./sum(d1(:));    

    dwarn=0;
    
    % --> if volume is too small, dilate until it is >= in size
    if sum(msk(:))<numelbnd(1)
        msk_new = msk;
        ndil=0;
        while sum(msk_new(:))<numelbnd(1) && dwarn==0
            %disp('mask too small! dilating!');
            ndil=ndil+1;
            if ndil>maxdil
                dwarn=1;
                warning('exceeded max dilation');
            else
                % box dilation of mask, constrained to admissible voxels...
                msk_new = double((convn( msk_new, ones(3,3,3), 'same') .* constr)>0);
            end        
        end
        % update mask object
        msk = msk_new; clear msk_new ndil; 
    end
    
    % --> if volume is too large, prune off enough voxels to ensure ==in size
    if sum(msk(:))>numelbnd(end)
        % distance from CoM per voxel
        [c,r,h] = meshgrid(1:ny,1:nx,1:nz);
        drch = sqrt( (r-mcom(1)).^2 + (c-mcom(2)).^2 + (h-mcom(3)).^2 );
        % get an ordered list of distances from CoM, and corresponding array indices
        di_inmsk(:,1) = drch(msk>0);
        di_inmsk(:,2) = sub2ind( [nx,ny,nz],r(msk>0),c(msk>0),h(msk>0) );
        % sort in increasing order, then keep up to numelbnd voxels
        di_inmsk = sortrows(di_inmsk,1); 
        di_inmsk = di_inmsk( 1:numelbnd(end), : );
        msk_new  = zeros(size(msk));
        msk_new(di_inmsk(:,2)) = 1;
        
        msk = msk_new; clear msk_new c r h drch di_inmsk;
    end
else
    dwarn = NaN;
end
% toc,

if sum(msk(constr==0))>0
    error('something went horribly wrong! resized mask is outside of constrained areas');
else
    finl_numvox = sum(msk(:));
end

%% (2) trim the vol to smallest cube that retains ROI > reduce unnecessary compute time 

% now get nonzero range of mask
ixx = find( squeeze(sum(sum(msk,2),3)) > 0 );
ixy = find( squeeze(sum(sum(msk,1),3)) > 0 );
ixz = find( squeeze(sum(sum(msk,1),2)) > 0 );
% trim mask and vol to fit range
vbox = vol(ixx(1):ixx(end),ixy(1):ixy(end),ixz(1):ixz(end));
mbox = msk(ixx(1):ixx(end),ixy(1):ixy(end),ixz(1):ixz(end));
% convert voxels outside roi to NaN, for compatibility with glca
vbox(mbox==0)=NaN;

% get texture - auto-set box limits, 
out = glca_kde( vbox, range_method, step_size, model, b_width, makefigs, 'default' );
% information on resizing -- catch if (a) start w wildly different sizes or (b) final size failed to hit bounds 
if isempty(numelbnd) 
    numelbnd = [NaN NaN];
elseif numel(numelbnd)==1
    numelbnd = [numelbnd numelbnd];
end
% [upper & low bounds][init/final/pct-cht][warn-]
out.sizecheck = [numelbnd, init_numvox, finl_numvox, 100*(finl_numvox-init_numvox)/init_numvox, dwarn];
