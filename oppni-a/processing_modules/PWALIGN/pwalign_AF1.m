function pwalign_AF1( Perffile, Basefile, prefix, odir1 )

% WARNING:>>> use odir2 to migrate necessary data for FWARP later?

% odir => sprintf('%s/warp',opath2f)
% prefix_set => sprintf('func%u',nr) x N_func
% Perffile => sprintf('%s/prewarp/func%u_2std.nii',opath2f,nr)  x N_func
% base_set => motref_0rel x N_func

% NB: when specifying alignment matrix outputs, make sure suffix is aff12.1D
%     throughout, otherwise cat_matvec will only do 1 line!

if ~exist( sprintf('%s/%s_aligned.nii.gz',odir1,prefix) ,'file')

    %-- a. rigid within-run alignment of all func volumes to func run-1 refbrick ( 
    unix(sprintf('3dvolreg -zpad 1 -base %s -1Dfile %s/%s_mpe -prefix %s/%s_pwalign.nii.gz -cubic %s', ...
        Basefile,odir1,prefix,odir1,prefix,Perffile));
end
