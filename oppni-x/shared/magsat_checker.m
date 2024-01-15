function magsat_checker( filelist )
% *
% . quick 'n dirty script to check whether initial scans are at
% . steady-state 
%

if ~isempty(dir('__opptmp_magsat_func*'))
    error('cannot run magsat - tempfiles suggest already in progress and/or early termination');
end

n_files = numel(filelist);

for n = 1:n_files

    [n, n_files],

    if ~exist(filelist{n},'file') 
        error('Cannot locate file: %s for magsat testing\n',filelist{n})
    end

    if contains( filelist{n}, '.nii.gz' )
        unix(sprintf('cp %s __opptmp_magsat_func.nii.gz',filelist{n}));
        unix('gunzip __opptmp_magsat_func.nii.gz');
    elseif contains( filelist{n}, '.nii' )
        unix(sprintf('cp %s __opptmp_magsat_func.nii',filelist{n}))
    else
        error('Unrecognized datatype of file: %s for magsat testing\n',filelist{n})
    end

    V = load_untouch_niiz('__opptmp_magsat_func.nii');
    volimg = double(V.img);
    volavg = mean(volimg,4);
    volimg = bsxfun(@rdivide, volimg, volavg);
    disp('!');
    masq = double( volavg > prctile(volavg(:),80));
    for t=1:size(volimg,4)
        %
        vtmp = volimg(:,:,:,t);
        tser(t,n) = mean(vtmp(masq>0));
    end

    % tempfile cleanup
    unix('rm __opptmp_magsat_func*')
end

% use fraction of subj as index -- more stable than mean change etc.
for(bsr=1:5000)
    list = ceil( n_files*rand(n_files,1) );
    avgV(:,bsr) = mean(  tser(1:30,list) ,2  );
    difF(:,bsr) = mean(  (tser(1:30,list)-tser(2:31,list))>0,2  );
    difV(:,bsr) = mean(  (tser(1:30,list)-tser(2:31,list)),2  );
end

figure, 
subplot(2,2,1);  
    plot( 1:30, mean(avgV,2), 'o-k', 'markerfacecolor',[0.5 0.5 0.5],'linewidth',2 ); hold on;
    plot( 1:30, prctile(avgV,[2.5 97.5],2), '.-k' ); hold on;
    plot([0 31],[1.0 1.0],'-r'); 
    ylim([0.99 1.01]); xlim([0 31]); title('mean intensity');
    ylabel('mean signal (standardized)'); xlabel('time (TR)')
subplot(2,2,3); 
    plot( 1:30, mean(difF,2), 'o-k', 'markerfacecolor',[0.5 0.5 0.5],'linewidth',2 ); hold on;
    plot( 1:30, prctile(difF,[2.5 97.5],2), '.-k' ); hold on;
    plot([0 31],[0.5 0.5],'-r'); 
    xlim([0 31]); title('fraction changed');
    ylabel('frac diff'); xlabel('time (TR)')
subplot(2,2,4);  
    plot( 1:30, mean(difV,2), 'o-k', 'markerfacecolor',[0.5 0.5 0.5],'linewidth',2 ); hold on;
    plot( 1:30, prctile(difV,[2.5 97.5],2), '.-k' ); hold on;
    plot([0 31],[0 0],'-r'); 
    xlim([0 31]); title('mean change');
    ylabel('mean diff'); xlabel('time (TR)')
