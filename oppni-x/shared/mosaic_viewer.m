function mosaic_viewer( volname, ax, title, savename, colormapp, sym )

if nargin<5
    colormapp = 'bone';
end
if nargin<6
    sym=0;
end

if ischar(volname)
    V = load_untouch_niiz(volname);
    vol = double(V.img); %clear V;
else
    vol = volname;
end

if ax==1 % sagittal slices (along x)
    vol = flip( permute( vol, [3 2 1]), 1);
elseif ax==2 % coronal slices (along y)
    vol = flip( permute( vol, [3 1 2]), 1);
else % axial slices (along z)
    %no flips needed
    vol = flip( permute( vol, [2 1 3]), 1);
end

br = squeeze(mean(mean(abs(vol),1),2));
istart = find( br>prctile(br,10),1,'first');
iend   = find( br>prctile(br,10),1,'last');

if (iend-istart+1) > 44
    ilist = unique( round(linspace(istart,iend,44)) );
else
    ilist = istart:iend;
end

f = figure;
f.Position=[100 475 2100 800]; % left edge, top edge, width, height
mosaic_big = [];
for q=1:4
    if numel(ilist)>((q-1)*11) 
        mosaic=[];
        for i=((q-1)*11+1):q*11
            if i<=numel(ilist)
            mosaic = [mosaic, vol(:,:,ilist(i))];
            else
            mosaic = [mosaic, zeros(size(vol,1),size(vol,2))];
            end
        end
        mosaic_big = [mosaic_big; mosaic];
    end
end
if sym==0
    bnd = prctile( mosaic_big(mosaic_big>eps), [1 99] );
else
    bnd = prctile( abs(mosaic_big(abs(mosaic_big)>eps)), 99 );
    bnd = [-bnd bnd];
end
imagesc( mosaic_big, bnd ); 
colormap(colormapp); colorbar;
set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
