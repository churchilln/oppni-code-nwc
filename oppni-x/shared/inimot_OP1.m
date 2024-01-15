function motref_0rel = inimot_OP1( Funcfile, prefix, odir )

% odir=> opath1f/init_mot_estim
% prefix => func(nr)
%
% --> checks for consistency of inimot!

% prefix for temp files
pref = [odir,'/__opptmp_p2func_inimot'];

if ~exist(sprintf('%s/%s_motref_0rel.txt',odir,prefix),'file')

    % build directory struct recursively
    unix(sprintf('mkdir -p %s',pref));

    % preparing data for first estimates of displacement...smooth 'n' mask
    unix(sprintf('3dmerge -prefix %s/func_smo6.nii -doall -1blur_fwhm 6 %s',pref,Funcfile));
    unix(sprintf('3dAutomask -prefix %s/func_smo6_mask.nii %s/func_smo6.nii',pref,pref));

    %--- a. mindisp estimation for alignment... (we can keep this, since it probably won't change after slice-proc)

    VF     = load_untouch_niiz( sprintf('%s/func_smo6.nii',pref) );
    MF     = load_untouch_niiz( sprintf('%s/func_smo6_mask.nii',pref) );
    % convert to 4D fMRI data volume (VV.img) into 2D matrix (voxels x time) for analysis
    epimat = nifti_to_mat( VF,MF ); 
    epimat = bsxfun(@minus,epimat,mean(epimat,2)); % mean-centering
    [u,l,v]=svd( epimat,'econ' );
    Qdat = (v(:,1:end-1)*l(1:end-1,1:end-1))';
    Qmed = median( Qdat,2 );
    Dist = sqrt(sum(bsxfun(@minus,Qdat,Qmed).^2));
    for t1=1:size(Qdat,2)
        for t2=1:size(Qdat,2)
            D2_rms(t1,t2) = sqrt( mean((Qdat(:,t1)-Qdat(:,t2)).^2) );
        end
    end
    % special cost-function vector (optimal is minimum)
    dmat(:,1) = (Dist-min(Dist))./(max(Dist)-min(Dist)); % standardized distance from medioid [0,1]
    dmat(:,2) = double( ([1; dmat(1:end-1,1)] > 0.5) | ([1;1; dmat(1:end-2,1)] > 0.5) | ([dmat(2:end,1) ;1] > 0.5) | ([dmat(3:end,1) ;1;1] > 0.5) ); % +1 if possible neighbouring spikes +- 2TR
    dmat(:,3) = zeros(numel(Dist),1); dmat(1:3,3)=1; dmat(end-2:end,3)=1; % +1 if start or end of the run
    [~,motref_1rel] = min( sum(dmat,2) ); % minimizer of 1, subject to costs 2 & 3
    motref_0rel = motref_1rel-1;

%     if exist( sprintf('%s/%s_motref_0rel.txt',odir,prefix),'file' )
%         xold = load( sprintf('%s/%s_motref_0rel.txt',odir,prefix) );
%         if xold ~= motref_0rel
%             error('motref seems to have changed - either restore param file, or delete P1 and try again with new setting!')
%         end
%     else
        % saving files
        dlmwrite( sprintf('%s/%s_motref_0rel.txt',odir,prefix), [motref_0rel] ); % ** found min-disp brick (zero-relative indexing)
        save(sprintf('%s/%s_dispdat.mat',odir,prefix),'dmat','motref_1rel','D2_rms');
%     end

    unix(sprintf('rm -rf %s',pref));

else
    disp('already compiled inimot')
    motref_0rel = load( sprintf('%s/%s_motref_0rel.txt',odir,prefix) );
end
