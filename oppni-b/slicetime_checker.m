function slicetime_checker( filelist )
%
% quick 'n dirty script to check slice tiing correction fideilty
%
% DEPRECATED .... DOESNT SEEM TO REALLY BE USEFUL!!!
%

n_files = numel(filelist);

for n = 1:n_files

    [n, n_files],

    if ~exist(filelist{n},'file') 
        error('Cannot locate file: %s for slicetime testing\n',filelist{n})
    end

    % takes zipped or unzipped -- defaults to gzip
    if contains( filelist{n}, '.nii.gz' )
        unix(sprintf('cp %s __opptmp_slicecheck_func.nii.gz',filelist{n}));
        unix('gunzip __opptmp_slicecheck_func.nii.gz')
    elseif contains( filelist{n}, '.nii' )
        unix(sprintf('cp %s __opptmp_slicecheck_func.nii',filelist{n}));
    else
        error('Unrecognized datatype of file: %s for slicechecker\n',filelist{n})
    end

    unix(sprintf('3dTshift -prefix __opptmp_slicecheck_func_tcr1.nii -tpattern alt+z __opptmp_slicecheck_func.nii'));
    unix(sprintf('3dTshift -prefix __opptmp_slicecheck_func_tcr2.nii -tpattern alt+z2 __opptmp_slicecheck_func.nii'));
    unix(sprintf('3dTshift -prefix __opptmp_slicecheck_func_tcr3.nii -tpattern alt-z __opptmp_slicecheck_func.nii'));
    unix(sprintf('3dTshift -prefix __opptmp_slicecheck_func_tcr4.nii -tpattern alt-z2 __opptmp_slicecheck_func.nii'));

    V = load_untouch_niiz('__opptmp_slicecheck_func.nii');
    volimg = double(V.img);
    volavg = mean(volimg,4);
    masq = double( volavg > prctile(volavg(:),80));
    mergmasq = masq(:,:,1:end-1) .* masq(:,:,2:end);
    zix = find( squeeze( sum(sum(mergmasq,1),2) ) > 20 );
    %
    volimg = volimg-mean(volimg,4);
    volimg = volimg./sqrt(sum(volimg.^2,4));

    for i=1:4
        V = load_untouch_niiz(sprintf('__opptmp_slicecheck_func_tcr%u.nii',i));
        volimg_t = double(V.img);
        volimg_t = volimg_t-mean(volimg_t,4);
        volimg_t = volimg_t./sqrt(sum(volimg_t.^2,4));
        %
        corset=[];
        for z=1:numel(zix)
            corm = sum( volimg_t(:,:,zix(z),:) .* volimg_t(:,:,zix(z)+1,:), 4 );
            corset = [corset; corm(mergmasq(:,:,zix(z))>0)];
        end
        cormed(n,i) = median( abs(corset) );
    end

    unix('rm __opptmp_slicecheck_func*')
end

figure, imagesc( cormed );
