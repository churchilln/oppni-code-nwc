function counts = clust_up_counter( vol )
%
vol = squeeze(vol);
if dim(vol,4)>1
    error('too many dims');
end

D=size(vol);              %% volume dimensions
% binarize the volume for cluster analysis
bino = double( abs(vol)>eps );

C1   = knn( bino );
C1   = C1.*bino;
C2   = knn(C1);
C2   = C2.*bino;
% at least 2 neighbours
bino_th   = double(C2>=2);

if( minclust <3 )
    error('minimum cluster size threshold of 3!');
elseif( sum(bino_th(:))==0 )
    disp('no clusters...returning empty');
    vol2 = zeros(size(vol));
else
    % now fit to minclust
    vox_idx      = find( bino_th > 0 );
    clust_set{1} = [vox_idx(1)];
    Nclust=1;

    for(i=2:length(vox_idx))
       %[i length(vox_idx)],    

       kn=zeros(Nclust,1);
       for(j=1:Nclust)
           [ic jc kc] =ind2sub(D, clust_set{j});
           [ix jx kx] =ind2sub(D, vox_idx(i));
           kn(j,1) = max( (abs(ic-ix)<=1) .* (abs(jc-jx)<=1) .* (abs(kc-kx)<=1) );
       end
       assn = find( kn );

       if( isempty(assn) )
           Nclust=Nclust+1;
           clust_set{Nclust} = vox_idx(i);
       elseif( length(assn)==1 )
           clust_set{assn} = [clust_set{assn}; vox_idx(i)]; 
       else
           % merging
           Nclust = Nclust - length(assn)+1;
           new_clust=[];
           for(k=1:length(assn)) new_clust = [new_clust; clust_set{assn(k)}]; end
           clust_set(assn)=[];
           clust_set{Nclust} = [new_clust; vox_idx(i)];
       end
    end

    for(i=1:Nclust) 
        ClustSize(i,1) = length(clust_set{i});
    end
end

counts = ClustSize;

%%
function count = knn( X )

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
