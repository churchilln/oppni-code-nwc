function out = fmri_interpolator( volname, outlier_struct, outlier_type, outprefix )
%
% now interpolates
%
%

% load fmri data
VV = load_untouch_niiz(volname);
VV.img = double(VV.img); %% format as double, for compatibility

if strcmpi(outlier_type,'0') || strcmpi(outlier_type,'OFF')
    outlier_vect = 0;
elseif strcmpi(outlier_type,'M')
    outlier_vect = outlier_struct.outl_mot;
elseif strcmpi(outlier_type,'V')
    outlier_vect = outlier_struct.outl_vol;
elseif strcmpi(outlier_type,'S')
    outlier_vect = outlier_struct.outl_slc;
elseif strcmpi(outlier_type,'MV') || strcmpi(outlier_type,'VM') || strcmpi(outlier_type,'1') || strcmpi(outlier_type,'ON')
    outlier_vect = outlier_struct.outl_vol_mot;
elseif strcmpi(outlier_type,'MS') || strcmpi(outlier_type,'SM')
    outlier_vect = outlier_struct.outl_slc_mot;
else
    error('unrecognized outlier choice for interpolation.')
end

if( sum(outlier_vect(:))==0 )
    disp('No outliers removed.');
else
    
    % dimensions of imaging data
    [Nx,Ny,Nz,Nt]=size(VV.img);
    TimeList = 1:Nt;
    
    % formatting the censor array
    sz = size(outlier_vect); % needs to be NtxNz
    if(min(sz)==1) 
        disp('single vector - applying same to all slices');
        X_cens = repmat(outlier_vect(:),[1 Nz]);
    elseif( sz(1)==Nt && sz(2)==Nz )
        disp('matrix - applying different adjustment to each slice');
        X_cens = outlier_vect;
    elseif( sz(1)==Nz && sz(2)==Nt )
        disp('matrix transposed - flipping then applying different adjustment to each slice');
        X_cens = outlier_vect';
    else
        error('problem with censor array?');
    end
    
    for(z=1:Nz)
        
        % slice-specific censor vec
        X_slcCens = X_cens(:,z); 
        
        % if outliers found, we drop + interpolate
        if( ~isempty(find(X_slcCens==1,1,'first')) )
            
            idx  = find(X_slcCens==0); % index of uncensored points
            
            % if point 1=outlier, make it same as first non-outlier scan
            if( X_slcCens(1) == 1 ) 
                X_slcCens(1)  = 0;
                VV.img(:,:,z,1)  = VV.img(:,:,z,idx(1)  );
            end
            % if endpoint=outlier, make it same as last non-outlier
            if( X_slcCens(Nt)== 1 )
                X_slcCens(Nt) = 0;
                VV.img(:,:,z,Nt) = VV.img(:,:,z,idx(end));
            end

            % take slice from 4D fMRI volumes; convert to vox x time matrix
            slcmat = reshape( VV.img(:,:,z,:), [],Nt );
            % get list of non-outlier timepoints, and image vectors
            subTimeList = TimeList(X_slcCens==0);
            submat     = slcmat(:,X_slcCens==0)';
            % interpolation -> reconstruct the full data matrix for all timepoints in (TimeList), 
            % based on the non-outlier input data (subTimeList, subvolmat)
            interpmat = interp1( subTimeList, submat, TimeList,'pchip' );
            % flip the volumes, so that it is (voxels x time)
            interpmat = interpmat';
            % plug back into matrix
            VV.img(:,:,z,:) = reshape( interpmat, Nx,Ny,1,Nt );            

        end
    end
end

% save the output results
save_untouch_niiz(VV,outprefix);   





