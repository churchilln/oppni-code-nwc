function amask_HY1(Adataset, Basedset, odir, ParamCell)
%
% .amask_HY1:
% .anatomical masking using hybrid methods

epath = mfilename('fullpath');
[ea,eb] = fileparts(epath);

if isempty(ParamCell)
    error('hybrid mask needs a list of submasks');
else

    if ~exist( sprintf('%s/anatBrainMask.nii.gz',odir) ,'file')

        nmask = numel(ParamCell);
    
        el_list = [];
        for n=1:nmask
            % path to next subfunc
            enam_sub = sprintf('amask_%s',ParamCell{n});
        
            % get function handle for analysis model of interest
            currPath=pwd;                               % get current path
            cd( ea );               % jump to module directory
            pfun= str2func( enam_sub ); % get function handle
            cd(currPath);                               % jump back to current path
            % execute step:
            pfun( Adataset, Basedset, odir, {} );  
            
            % now rename copy of file
            unix(sprintf('mv %s/anatBrainMask.nii.gz %s/anatBrainMask_meth%u.nii.gz',odir,odir,n))
            el_list = [el_list, ' ', sprintf('%s/anatBrainMask_meth%u.nii.gz',odir)];
        end
    
        % now get the intersection
        masksum = 0;
        for n=1:nmask
            M = load_untouch_niiz(sprintf('%s/anatBrainMask_meth%u.nii.gz',odir,n));
            masksum = masksum + double(M.img);
        end
        M.img = double(masksum==nmask);
        save_untouch_niiz(M,sprintf('%s/anatBrainMask.nii.gz',odir));
        
    else
        disp('hybrid-mask already exists!')
    end
end
