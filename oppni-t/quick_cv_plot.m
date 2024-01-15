close all;
clear;

e = dir('texturial_matfiles_optim/sub-*.mat');

% for j=1:13
% figure;
% for i=1:16
% 
%     load([e(i).folder,'/',e(i).name]);
%     subplot(4,4,i); plot( log(out.SZfull(:,:,j)), zscore(out.CVfull(:,:,j)), '.-' );
% 
% end
% end

for i=1:numel(e)

    load([e(i).folder,'/',e(i).name]);
    for j=1:13 % indexing>direction
        for k=1:160 % indexing>roi
            [vx,ix] = min( out.CVfull(:,j,k) );
            szarr(i,j,k) = out.SZfull(ix,j,k);
            cvarr(i,j,k) = vx;
            ixarr(i,j,k) = ix;
            nvarr(i,k) = out.sizecheck(k,6);
        end
    end
end

bnd = prctile( szarr(:), [5 95]);
for i=1:12
    blk = permute( szarr(:,i,:), [1 3 2]);
    subplot(3,4,i); imagesc( blk, bnd )
end

return;


figure;
for k=1:13
    subplot(2,7,k); imagesc( szarr(:,:,k), prctile(szarr(:),[2.5 97.5]) ); colormap jet;
end

figure;
for k=1:13
    subplot(2,7,k); imagesc( cvarr(:,:,k), prctile(cvarr(:),[2.5 97.5]) ); colormap jet;
end
