function make_effex_maps( file0, file1, mask, type )
%
% .this script generates a summative map of change from file0 to file1
% .includes: varrat, connloc, connglob, 
%

V0 = load_untouch_niiz(file0);
V1 = load_untouch_niiz(file1);
M  = load_untouch_niiz(mask);

m = double(M.img);

vmat0 = nifti_to_mat(V0,M);
vmat1 = nifti_to_mat(V1,M);

if strcmpi(type,'varrat')

    img = var( double(V1.img), 0,4) ./ (var( double(V0.img), 0,4) + eps);
    img = img .* double(M.img);
    
elseif strcmpi(type,'connloc')

    a = double(V0.img) .* m;
    a = bsxfun(@minus,a,mean(a,4)); % mean center
    a = bsxfun(@rdivide,a,sqrt(sum(a.^2,4))); % unit norm

    % sum correlations with all spatially shifted versions
    img = 0;
    img = img + sum( a .* circshift(a, 1,1), 4);
    img = img + sum( a .* circshift(a,-1,1), 4);
    img = img + sum( a .* circshift(a, 1,2), 4);
    img = img + sum( a .* circshift(a,-1,2), 4);
    img = img + sum( a .* circshift(a, 1,3), 4);
    img = img + sum( a .* circshift(a,-1,3), 4);
    % sum number of times this gives in-mask entry
    msq = 0;
    msq = msq + ( m .* circshift(m, 1,1) );
    msq = msq + ( m .* circshift(m,-1,1) );
    msq = msq + ( m .* circshift(m, 1,2) );
    msq = msq + ( m .* circshift(m,-1,2) );
    msq = msq + ( m .* circshift(m, 1,3) );
    msq = msq + ( m .* circshift(m,-1,3) );

    img = img ./ (msq + eps);
end
