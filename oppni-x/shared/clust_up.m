function vol2 = clust_up( vol, minclust, vol_select, into_parc,out_format, file_name,coord_spc )
%
% =========================================================================
% CLUST_UP: Cluster-size thresholding on 3D volume
% =========================================================================
%
% Syntax:
%           vol2 = clust_up( vol, minclust, (vol_select,into_parc,out_format,file_name,coord_spc) )
%
% Input:
%           vol: 3D volume, clusterized based on binary thresholding
%                (x==0 vs x~=0)
%      minclust: minimum cluster-size threshold. Must be 3 or more
%
%    vol_select: integer (any value >=1). index selecting volume, in case of 4D volume set (default=1)
%     into_parc: binary value (0,1). export cluster-corrected image as set of parcel labels? (default=0)
%    out_format: string ('none','standard','simple'). specifies format of output file (default='none')
%     file_name: string. name of file for output report (default='clust_up_output.txt')
%     coord_spc: string ('A','V'). coordinates of output. Either  as array index values ('A') or voxel-space coordinates ('V'). (default='A')
%
% Output:
%          vol2: cluster thresholded image (original or parcel labelled)
%

if nargin<2
    error('need to at least specify volume and minimum cluster size')
end
if nargin<3 || isempty(vol_select) % which volume (if 4D series)
    vol_select = 1;
end
if nargin<4 || isempty(into_parc) % break into indexed parcels?
    into_parc = 0;
end
if nargin<5 || isempty(out_format) % break into indexed parcels?
    out_format = 'none';
end
if nargin<6 || isempty(file_name) % break into indexed parcels?
    file_name = 'clust_up_output.txt';
end
if nargin<7 || isempty(coord_spc) % coordinates for output
    coord_spc = 'A';
end

if(ischar(vol))
    % load n' clear
    volnii = vol; % save filename
    V=load_untouch_niiz(vol);
    vol = double(V.img);
    fromnii = true;
    if strcmpi(coord_spc,'V')
    cubevol = prod( V.hdr.dime.pixdim(2:4) );
    units   = 'mm^3';
    elseif strcmpi(coord_spc,'A')
    cubevol = 1;
    units   = 'vox';
    else
    error('no recognizze coord_spc')
    end
else
    fromnii = false;
    cubevol = 1;
    units   = 'vox';
    %--
    coord_spc='A'; % force it to take array indices
end
if strcmpi(out_format,'none')
    toout=0;
else
    toout=1;
end
if strcmpi(coord_spc,'V')
    omat(1,:) = V.hdr.hist.srow_x;
    omat(2,:) = V.hdr.hist.srow_y;
    omat(3,:) = V.hdr.hist.srow_z;
    om1 = omat(1:3,1:3);
    if sum(sum(( om1 - diag(diag(om1)) )))>eps
        error('non-cardinal axes in header - not currently supported!')
    else
        ov_scal = diag(om1);
        clear om1;
    end
    ov_trns = omat(1:3,4);
    clear omat;
end

vol = squeeze(vol);  %% squash out singleton dimensions
if size(vol,5)>1
    error('input data structure is too complex! trim it down!');
end
D=size(vol);  %% volume dimensions

if(numel(D)>3 && D(4)>1 ) %% if a 4D timeseries
    if    ( isempty(vol_select) ) error('4D series. need to select a volume');
    elseif( vol_select>D(4)     ) error('volume selection index exceeds number of volumes');
    elseif( numel(vol_select)>1 ) error('can only select a single volume for cluster analysis');
    else
        % select single volume, adjust dimensions
        vol=vol(:,:,:,vol_select);
        D(4)=1;
    end
elseif( length(D)<3) 
    D(3)=1;
end

% binarize the volume for cluster analysis
bino = double( abs(vol)>eps );
% count neighbours + neighbours-of-neighbours
CT   = knn( bino, 1 );
CT   = knn(CT, 1);
% binarize to keep only if at least 2 neighbours
bino   = double(CT>=2);

if( minclust <3 )
    error('minimum cluster size threshold of 3!');
elseif( sum(bino(:))==0 )
    disp('no clusters...returning empty');
    vol2 = zeros(size(vol));
else
    % now fit to minclust

    vox_idx      = find( bino > 0 ); % find all nonzero voxels
    clust_set{1} = [vox_idx(1)]; % initialize with one cluster, containing first nonzero vox 
    Nclust=1;

    for(i=2:length(vox_idx)) % pick through remaining nonzero vox

       [ix jx kx] =ind2sub(D, vox_idx(i)); % array indices

       kn=zeros(Nclust,1); % to thru ea cluster, find if it's connected to any vox
       for(j=1:Nclust)
           [ic jc kc] =ind2sub(D, clust_set{j}); % array indices of all target vox
           kn(j,1) = max( (abs(ic-ix)<=1) .* (abs(jc-jx)<=1) .* (abs(kc-kx)<=1) );
       end
       assn = find( kn ); % find nonzero entries

       if( isempty(assn) ) % if unconnected to existing, append a new cluster
           Nclust=Nclust+1;
           clust_set{Nclust} = vox_idx(i);
       elseif( length(assn)==1 ) % if connected to 1 existing, add to cluster
           clust_set{assn} = [clust_set{assn}; vox_idx(i)]; 
       else % if connected to multiple clusters...merge all of them!
           % reduce number of clusters
           Nclust = Nclust - length(assn)+1; 
           % put all voxels in the connected clusters into a new one
           new_clust=[];
           for(k=1:length(assn)) new_clust = [new_clust; clust_set{assn(k)}]; end
           % now delete the old subclusters
           clust_set(assn)=[];
           % put the new voxel into the merged cluster
           clust_set{Nclust} = [new_clust; vox_idx(i)];
       end
    end
    disp('done populating clusters!');

    % cluster sizes
    for(i=1:Nclust) 
        ClustSize(i,1) = length(clust_set{i});
    end
    % reorder clusters by decreasing size
    rowix = sortrows([(1:Nclust)' ClustSize],-2);
    clust_set=clust_set(rowix(:,1));
    % now delete below-threshold clusters
    clust_map=zeros(size(vol));
    Nclust2=0;
    for(i=1:Nclust)
        if( length(clust_set{i}) >= minclust )
            Nclust2 = Nclust2+1; 
            clust_map(clust_set{i})=Nclust2;
        end
    end
    
    if(toout)
        % initialize output file
        fin = fopen(file_name,'wt');
        fprintf(fin,'\nCluster-size report, minimum cluster size =%s',num2str(minclust));
        if(fromnii) fprintf(fin,'\ninput files =%s\n\n',volnii);
        else        fprintf(fin,'\nused 3d input matrix\n\n');
        end
    end
    
    strout = ['total number of clusters: ',num2str(Nclust2)];
    disp(strout);
    if(toout) fprintf(fin,[strout,'\n\n']); end

    % get some relevant stats on clusters
    for(n=1:Nclust2)

        vtemp  = clust_map(clust_map(:) == n);
        bintmp = double(clust_map==n);
        voltmp = vol.*double(clust_map==n);
        %%%
        d1=sum(sum(abs(voltmp),2),3); xc = sum(d1(:).*(1:size(bintmp,1))')./sum(d1(:));
        d1=sum(sum(abs(voltmp),1),3); yc = sum(d1(:).*(1:size(bintmp,2))')./sum(d1(:));
        d1=sum(sum(abs(voltmp),1),2); zc = sum(d1(:).*(1:size(bintmp,3))')./sum(d1(:));    
        val{1}=round([xc yc zc]); % center of mass
        val{2} = numel(vtemp); %cluster size
        %%%
        [x y z] = ind2sub( D, find(abs(voltmp)==max(abs(voltmp(:)))));
        val{3} = [x,y,z]; %peak coordinates
        val{4} = voltmp(x(1),y(1),z(1)); %peak value
        
        % adjust coordinates?
        if( strcmpi(coord_spc,'V') )
            % 
            val{1} = ((val{1} - 1) .* ov_scal(:)') + ov_trns(:)';
            val{3} = ((val{3} - 1) .* ov_scal(:)') + ov_trns(:)';
        end
        
        strout = ['Clust#',num2str(n),...
                   ': size(',units,')=',num2str(cubevol*val{2}),...
                   ', CoM (',coord_spc,') =[',num2str(round(val{1}(1))),',',num2str(round(val{1}(2))),',',num2str(round(val{1}(3))),']',...
                   ', peak value  =',sprintf('%0.2f',val{4}),...
                   ', peak (',coord_spc,') =[',num2str(round(val{3}(1))),',',num2str(round(val{3}(2))),',',num2str(round(val{3}(3))),']',...
                   ];
        disp(strout);
        if(toout)
            if( strcmpi(out_format,'standard') )
                fprintf(fin,[strout,'\n']); 
            elseif(  strcmpi(out_format,'simple') )
                fprintf(fin,'%u, %u,%u,%u, %f, %u,%u,%u\n',cubevol*val{2}, val{1}(1),val{1}(2),val{1}(3), val{4}, val{3}(1),val{3}(2),val{3}(3));
            else
                error('type not recognized');
            end
        end
    end
    
    if(into_parc==0)
        vol2 = vol.*double(clust_map>0);
        parcstr='';
    else
        vol2 = clust_map;
        parcstr=['_Parcels',num2str(Nclust2)];
    end
    
    if(toout) fclose(fin); end
    
    if(fromnii)
        [~,pref,ext] = fileparts(volnii);
        V.img = vol2;
        V.hdr.dime.dim(5:end) = 1;
        if( into_parc<2 )
            save_untouch_niiz(V,[pref,'_clustSize',num2str(minclust),'',parcstr,'.nii']);
        else
            for(i=1:Nclust2)
                V.img = double(vol2==i);
                save_untouch_niiz(V,[pref,'_clustSize',num2str(minclust),'_P',num2str(i),'of',num2str(Nclust2),'.nii']);
            end
        end
    end
end
% % nii=M;
% % nii.img = tmp;
% % nii.hdr.dime.datatype = 16;
% % nii.hdr.hist = V.hdr.hist;
% % nii.hdr.dime.dim(1) = 4;
% % nii.hdr.dime.dim(5) = 1;
% % save_untouch_niiz(nii,['recov_results/',modlab{modal},'_lng2.nii']);
%%
function count = knn( X, occup_only )
% counter number of "occupied" neighbours at each voxel
% for binary array X

D=size(X); if(length(D)<3) D(3)=1; end
count = 0;
count = count + cat( 1, zeros([1,D(2),D(3)]), X(1:end-1,:,:) );
count = count + cat( 1, X(2:end,:,:), zeros([1,D(2),D(3)]) );
count = count + cat( 2, zeros([D(1),1,D(3)]), X(:,1:end-1,:) );
count = count + cat( 2, X(:,2:end,:), zeros([D(1),1,D(3)]) );
if(D(3)>1)
count = count + cat( 3, zeros([D(1),D(2),1]), X(:,:,1:end-1) );
count = count + cat( 3, X(:,:,2:end), zeros([D(1),D(2),1]) );
end

if occup_only>0
    count = count .* X;
end