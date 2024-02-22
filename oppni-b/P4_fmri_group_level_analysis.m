function P4_fmri_group_level_analysis( inputfile, pipelinefile, paramlist, outpath, param, volno, design_mat, analysis_model, model_contrast, THRESH_METHOD, out_folder_name, censor )
%
% Input:
%
%          inputfile : string, giving name of input textfile listing data to process
%      pipelinefiles : string, giving name of pipeline textfile specifying which pipeline steps to apply to the data
%          paramlist : string, giving name of parameter file specifying special arguments to use in pipelines
%            outpath : string, specifying destination directory for processed outputs
%              param : string, specifying which derived param to analyze 
%              volno : index number among contrasts to analyze, if multiple ones
%                      can also use this to specify contrast
%         design_mat : design matrix in table format; should have header w variable names
%     analysis_model : string specifying analysis model ('Ttest' or 'GLM')
%                      (default -- Ttest)
%     model_contrast : string specifying contrast, based on design matrix fields
%                     (default -- none, do a 1-sample contrast)
%      THRESH_METHOD : threshold method -- 2-element cell array specifying type, critical value}
%                      can be uncorrected: {'UNCORR',[p-value]}
%                      fdr-corrected: {'FDR',[q-value]}
%                      or cluster-size corrected: {'CLUST',[minimum cluster size]}
%                      (default -- none, dont create a thresholded map)
%    out_folder_name : string specifying prefix for output files
%                      format is ['Analysis_',<outname>]
%             censor : file tells you which participants to RETAIN (1=keep /0=exclude
%


% initializing structure and running checkes
[subject_list, ~, PipeStruct_aug, ParamStruct_aug] = P0_fmri_populateDirectories( inputfile, pipelinefile, paramlist, outpath );

if isempty(censor)
    censor = ones(numel(subject_list),1);
elseif ischar(censor)
    if exist(censor,'file')
        censor = readtable(censor);
        censor = table2array(censor(:,end));
    else
        error('cannot find censor file');
    end
else
    if size(censor,2)>1
        error('numeric input for censor must be a vector')
    end
end

% now augmenting outpath... do this after P0!
outpath = fullfile(outpath,'fmri_proc'); % subdir should be fmri_proc

%% Loading imaging data

MBstr = ([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_brain_mask_grp.nii']);
MB = load_untouch_niiz(MBstr);
maskS = double(MB.img);

for i=1:numel(subject_list)

    filepath = [outpath,'/',subject_list{i},'/func_proc_p2/Pipe_',PipeStruct_aug.PNAME{1},'/',ParamStruct_aug.Variable_ID,'/out_analysis.mat'];
    x=load(filepath);

    if i==1
        fi = fieldnames( x.out_analysis.image );
        if sum(strcmpi( param, fi))>0
            param_type = 'image';
            matdims = [size(x.out_analysis.image.(param),1) 1];
        else
            fi = fieldnames( x.out_analysis.mat2d );
            if sum(strcmpi( param, fi))>0
                param_type = 'mat2d';
                matdims = size(x.out_analysis.mat2d.(param));
            else
                error('unrecognized param type??')
            end
        end
    end

    if strcmpi( param_type,'image' )
        if isnumeric(volno)
            if numel(volno)==1
                datamat(:,i) = x.out_analysis.image.(param)(:,volno);
            else
                datamat(:,i) = x.out_analysis.image.(param) * volno(:);
            end
        else
            volno2 = regexp( volno,'-','split');
            bb = str2num( strrep(volno2{1},'+',',') );
            aa = str2num( strrep(volno2{2},'+',',') );
            datamat(:,i) = mean(x.out_analysis.image.(param)(:,bb),2) - mean(x.out_analysis.image.(param)(:,aa),2);
        end
    elseif strcmpi( param_type,'mat2d' )
        if ~isempty(volno) error('condition constrasts unsupported for connectivity arrays'); end
        datamat(:,i) = reshape( x.out_analysis.mat2d.(param), [],1);
    else
        error('??')
    end
end
% applying censoring to datamat
D2 = datamat(:,censor==1);

%% Loading design matrix / behavioural data

if isempty(design_mat)
    Xdes = [];
    heads = [];
else
    if ischar(design_mat)
        if exist(design_mat,'file')
            design_mat = readtable(design_mat);
            % trim off the first col., save as participant ID list
            pID_des    = design_mat(:,1);
            design_mat = design_mat(:,2:end);
        else
            error('cannot find design matrix file');
        end
    else
        error('requires table-format design matrix for now!');
    end
    disp('design matrix loaded!');
    if size(design_mat,1) ~= size(datamat,2)
        error('design matrix size does not match number of fMRI participants loaded!');
    end
    Xdes = table2array(design_mat);
    heads = design_mat.Properties.VariableNames;    
end

if ~isempty(model_contrast) % split into strings
    model_contrast = regexp( model_contrast,',','split');
    is_interact = zeros(numel(model_contrast),1);
    for i=1:numel(model_contrast)

        if contains(model_contrast{i},'*')
            model_contrast_sub = regexp( model_contrast{i},',','split');
            ix1 = find( strcmpi( strtrim(model_contrast_sub{1}), heads) );
            if isempty(ix1) error('cannot find variable "%s" in design matrix!\n',model_contrast_sub{1}); end
            ix2 = find( strcmpi( strtrim(model_contrast_sub{2}), heads) );
            if isempty(ix2) error('cannot find variable "%s" in design matrix!\n',model_contrast_sub{2}); end

            Xdes_new(:,i) = Xdes(:,ix1) .* Xdes(:,ix2);
            %
            Xdes_new_sub(:,i,1) = Xdes(:,ix1);
            Xdes_new_sub(:,i,2) = Xdes(:,ix2);
            is_interact(i) = 1;
        else
            ix = find( strcmpi( strtrim(model_contrast{i}), heads) );
            if isempty(ix) error('cannot find variable "%s" in design matrix!\n',model_contrast{i}); end
            Xdes_new(:,i) = Xdes(:,ix);
            %
            Xdes_new_sub(:,i,1:2) = zeros(size(Xdes,1),1,2);
            is_interact(i) = 0;
        end
    end
    % applying censoring to design matrix
    X2 = Xdes_new(censor==1,:);
    X2_sub = Xdes_new_sub(censor==1,:,:);
else
    Xdes_new = [];
    Xdes_new_sub = [];
    X2 = [];
    X2_sub = [];
end

%% Check #1 >> missing data?

fprintf('\n=========================================================\n');

    % check for bad/missing values
    fprintf('Missing data check, imaging:\n');
    ixdrop_D = (mean(~isfinite(D2),1)>0.10)';
    fprintf('number of volumes with more than 10 percent of voxels missing: %s\n',sum(ixdrop_D));
    if     mean(ixdrop_D)==0   fprintf('no missing data.')
    elseif mean(ixdrop_D)<0.10 fprintf('less than 10 percent of participants have substantial missing data');
    elseif mean(ixdrop_D)>0.10 warning('more than 10 percent of participants have substantial missing data');
    elseif mean(ixdrop_D)>0.20   error('more than 20 percent of participants have substantial missing data. Halting!');
    end
    fprintf('Missing data check, design matrix:\n');
    if isempty(X2)
        disp('...no design matrix!\n');
    else
        ixdrop_X = ~isfinite(mean(X2,2));
        fprintf('number of design rows with missing data: %s\n',sum(ixdrop_X));
        if     mean(ixdrop_X)==0   fprintf('no missing data.')
        elseif mean(ixdrop_X)<0.10 fprintf('less than 10 percent of participants have substantial missing data');
        elseif mean(ixdrop_X)>0.10 warning('more than 10 percent of participants have substantial missing data');
        elseif mean(ixdrop_X)>0.20   error('more than 20 percent of participants have substantial missing data. Halting!');
        end
    end
    fprintf('number of datapoints dropped, total: %u (from original %u)\n',sum(ixdrop_X | ixdrop_D),size(D2,2) );
    fprintf('dropping any missing points from data before proceeding...\n')
    %-- dropping  points with missing data before proceeding
    D2(:,ixdrop_X | ixdrop_D) = [];
    X2(ixdrop_X | ixdrop_D,:) = [];
    X2_sub(ixdrop_X | ixdrop_D,:,:) = [];

fprintf('=========================================================\n\n');

%% Check #2 >> quality of design matrix

fprintf('\n=========================================================\n');
    
    fprintf('Quality checks, design matrix:\n');
    if isempty(X2)
        fprintf('design matrix is empty, skipping this step....\n')
    else
        ix = find( sum(X2.^2,1) < eps );
        if ~isempty(ix)
            strtmp = []; for i=1:numel(ix) [strtmp, ', ', model_contrast{i}]; end
            error('empty columns in: %s\n',strtmp(3:end));
        else
            fprintf('no empty columns!\n');
        end
        xtmp = max(abs(X2),[],1);
        if numel(xtmp)>1 && max(xtmp)/min(xtmp) > 1E6
            [~,ix1]=min(xtmp);
            [~,ix2]=max(xtmp);
            error('Numeric scaling issue likely to give poor results: check %s or %s (possibly others)\n',model_contrast{i1},model_contrast{i2});
        else
            fprintf('scaling ok!\n')
        end
        % check indiv cols
        for i=1:size(X2,2)
    
            strtmp = ['col.',num2str(i),' (',model_contrast{i},') -- '];

            if numel(unique(X2(:,i)))==2
                strtmp = [strtmp 'binary: '];
                ncl1= sum( X2(:,i)==min(X2(:,i)) );
                ncl2= sum( X2(:,i)==max(X2(:,i)) );
                strtmp = [strtmp, 'sample split: ',num2str(ncl1),'/',num2str(ncl2),'.'];
                if    ( ncl1/(ncl1+ncl2) < 0.10 )   strtmp = [strtmp ' grp1 constitutes <10 pct. of your sample. Extremely unbalanced!\n'];
                elseif( ncl2/(ncl1+ncl2) < 0.10 )   strtmp = [strtmp ' grp2 constitutes <10 pct. of your sample. Extremely unbalanced!\n'];
                elseif( ncl1/(ncl1+ncl2) < 0.20 )   strtmp = [strtmp ' grp1 constitutes 10-20 pct. of your sample. Unbalanced - proceed with caution!\n'];
                elseif( ncl2/(ncl1+ncl2) < 0.20 )   strtmp = [strtmp ' grp2 constitutes 10-20 pct. of your sample. Unbalanced - proceed with caution!\n'];
                else  strtmp = [strtmp ' ok.\n'];
                end       
            else
                strtmp = [strtmp 'continuous: '];
                strtmp = [strtmp, sprintf('skew=%.02f, kurt=%.02f.',skewness(X2(:,i)),kurtosis(X2(:,i)))];
                sknul = skewness(randn(size(X2,1),5000),0,1);
                isfa = 0;
                if mean( abs(skewness(X2(:,i))) > abs(sknul) ) > 0.95
                    isfa = 1;
                    strtmp = [strtmp ' variable is significantly non-normal (skew test, p<.05). Proceed with caution!\n'];
                end
                sknul = kurtosis(randn(size(X2,1),5000),0,1);
                if mean( abs(kurtosis(X2(:,i))) > abs(sknul) ) > 0.95
                    isfa = 1;
                    strtmp = [strtmp ' variable is significantly non-normal (kurt test, p<.05). Proceed with caution!\n'];
                end 
                if isfa==0
                    strtmp = [strtmp ' ok.\n'];
                end
            end
            fprintf(strtmp);
        end
        % collinearity
        cctmp = corr(X2);
        cctmp = cctmp .* triu( ones(size(cctmp)), 1);
        [vx,ix] = max(abs(cctmp(:)));
        if numel(xtmp)>1 && vx>0.75
            [i1,i2]=ind2sub(size(cctmp),ix);
            error('Redundant predictors: check %s or %s (possibly others!)\n',model_contrast{i1},model_contrast{i2});
        else
            fprintf('pairwise collin ok!\n')
        end
        % condition
        Xtmp = [ones(size(X2,1),1) X2]; % include intercept!
        cnum = cond(Xtmp);
        fprintf('Condition check. Number of regressors (excluding constant): %u. Condition number: %.03f.\n',size(X2,2),cnum );
        if    ( cnum > 1000 ) fprintf(['GLM design matrix is probably multi-collinear! condition number>1000\n']);
        elseif( cnum > 100  ) fprintf(['GLM design matrix is approaching multi-collinearity! condition number>100\n']);
        else                  fprintf(['GLM design matrix is probably ok! condition number<100\n']);
        end
    end

fprintf('=========================================================\n\n');


%% Check #2 >> quality of datamatrix

fprintf('\n=========================================================\n');
    
    figure;
    dtmp = D2 - mean(D2,2);
    subplot(3,2,1); imagesc( dtmp ); title('mean-centered')
    dtmp = dtmp./std(dtmp,0,2);
    subplot(3,2,2); imagesc( dtmp ); title('var-normed')
    dtmp = dtmp.^2;
    subplot(3,2,2); imagesc( dtmp ); title('squared deviation')
    
    dif = bsxfun(@minus,D2,mean(D2,2,'omitnan')).^2; % deviation from mean map
    outl = mean(dif,1)';                       % mean deviation (averaging over voxels)
    outl = outl./max(outl);                    % renorming deviation
    figure;
    subplot(3,2,2); bar( outl ); ylim( [0 1.01]);
    PARMHAT = gamfit(outl(isfinite(outl)));    % gamma distribution fitting 
    Pgam = gamcdf(outl,PARMHAT(1),PARMHAT(2)); % probability on gammas
    subplot(3,2,4); bar( 1-Pgam(isfinite(outl)) ); ylim( [0 0.05]);
    [p th]=fdr(1-Pgam,'p',0.05,0);             % signifiant outliers FDR=0.05
    if sum(th)<=0
        fprintf('no significant outlier volumes!\n');
    else
        fprintf('Found %u significant outlier volumes:\n',sum(th>0));
        fith = find(th>0);
        strtmp = [];
        for i=1:numel(fith)
            strtmp = [strtmp ', ', subject_list{fith(i)}];
        end
        fprintf('   %s\n',strtmp(3:end))
    end

fprintf('=========================================================\n\n');


if strcmpi( analysis_model, 'GLM' )
    
    if isempty(model_contrast) error('need to specify a model contrast for GLMs!'); end
    if isempty(X2) error('need to load a design matrix for GLMs!'); end
    
    out_analysis = GLM_ph( D2, X2 );

elseif strcmpi( analysis_model, 'GLMboot' )

    if isempty(model_contrast) error('need to specify a model contrast for GLMs!'); end
    if isempty(X2) error('need to load a design matrix for GLMs!'); end
    
    out_analysis = GLMboot_ph( D2, X2 );

elseif strcmpi( analysis_model, 'Ttest' )
    
    if isempty(X2)
        disp('no contrast. defaulting to 1-sample t-test');
    else
        if size(X2,2)>1 error('can only have a single (categorical) predictor for T-testing'); end
        if numel(unique(X2)) ~= 2
            error('for t-test, predictor must be categorical');
        end
    end

    out_analysis = Ttest_ph( D2, X2 );

elseif strcmpi( analysis_model, 'Bootstrap' )
    
    if isempty(X2)
        disp('no contrast. defaulting to 1-sample t-test');
    else
        if size(X2,2)>1 error('can only have a single (categorical) predictor for T-testing'); end
        if numel(unique(X2)) ~= 2
            error('for t-test, predictor must be categorical');
        end
    end

    out_analysis = Bootstrap_ph( D2, X2 );

elseif strcmpi( analysis_model, 'Permute' )
    
    if isempty(X2)
        disp('no contrast. defaulting to 1-sample t-test');
    else
        if size(X2,2)>1 error('can only have a single (categorical) predictor for T-testing'); end
        if numel(unique(X2)) ~= 2
            error('for t-test, predictor must be categorical');
        end
    end

    out_analysis = Permute_ph( D2, X2 );

elseif strcmpi( analysis_model, 'PCA' )

    figure,imagesc( zscore(X2')',[-3.5 3.5]); colormap jet;

    [u,l,v] = svd( X2,'econ' );

    figure,
    subplot(2,2,1); bar( diag(l.^2)./trace(l.^2) ); title('design matrix')
    subplot(2,2,2); plot( u(:,1), u(:,2), 'ok', 'markerfacecolor',[0.5 0.5 0.5] )
    subplot(2,2,3); bar( u(:,1:2) );
    subplot(2,2,4); bar( v(:,1:2) );

    [u,l,v] = svd( D2,'econ' );

    figure,
    subplot(2,2,1); bar( diag(l.^2)./trace(l.^2) ); title('design matrix')
    subplot(2,2,2); plot( v(:,1), v(:,2), 'ok', 'markerfacecolor',[0.5 0.5 0.5] )
    subplot(2,2,3); bar( v(:,1:2) );
    subplot(2,2,3); bar( v(:,1:2) );
    subplot(2,2,4); bar( u(:,1:2) );
    axial_plot( u(:,1), maskS, 6, 1, 1 );
    axial_plot( u(:,2), maskS, 6, 1, 1 );

    disp('just plotting for now---');
    return;
end

mkdir_r(out_folder_name);
out_analysis.submask = x.out_analysis.submask;
save([out_folder_name,'/out_analysis.mat'],'out_analysis');

a = fieldnames(out_analysis);
pix = find( contains(a,'_p'));
if isempty(pix) error('no p-value thresholding available!')
elseif numel(pix)>1 error('too many p-value fields');
end
pfield = a{pix};
pix = find( strcmpi(a,pfield(1:end-2)));
if isempty(pix) error('no effect size matching p-value!')
elseif numel(pix)>1 error('too many effect sizes matching p-value?');
end
vfield = a{pix};

%%

tmaps = out_analysis.(vfield);
pmaps = out_analysis.(pfield);

if strcmpi( param_type,'image' )
    
    %--1
    TMPVOL = zeros( [size(out_analysis.submask), size(tmaps,2)] );
    for i=1:size(tmaps,2)
        tmp = out_analysis.submask;
        tmp(tmp>0) = tmaps(:,i);
        TMPVOL(:,:,:,i) = tmp;
    end
    nii=MB;
    nii.img = TMPVOL;
    nii.hdr.dime.datatype = 16;
    nii.hdr.hist = MB.hdr.hist;
    nii.hdr.dime.dim(5) = size(tmaps,2);
    save_untouch_niiz(nii,[out_folder_name,'/',vfield,'_unthresh.nii']); 
    
    %--2
    TMPVOL = zeros( [size(out_analysis.submask), size(tmaps,2)] );
    for i=1:size(tmaps,2)
        tmp = out_analysis.submask;
        tmp(tmp>0) = pmaps(:,i);
        TMPVOL(:,:,:,i) = tmp;
    end
    nii=MB;
    nii.img = TMPVOL;
    nii.hdr.dime.datatype = 16;
    nii.hdr.hist = MB.hdr.hist;
    nii.hdr.dime.dim(5) = size(tmaps,2);
    save_untouch_niiz(nii,[out_folder_name,'/',pfield,'_unthresh.nii']); 
    
    if strcmpi(THRESH_METHOD{1},'FDR')
        [~,th] = fdr( pmaps,'p',THRESH_METHOD{2},0 );
        numaps = tmaps .* th;
    elseif strcmpi(THRESH_METHOD{1},'CLUST')
        for i=1:size(tmaps,2)
           tmp = out_analysis.submask;
           tmp(tmp>0) = double(pmaps(:,i)<=0.005) .* tmaps(:,i);
           tmp = clust_up( tmp, THRESH_METHOD{2} );
           numaps(:,i) = tmp(out_analysis.submask>0);
        end
    elseif strcmpi(THRESH_METHOD{1},'UNCORR')
        numaps = tmaps .* double(pmaps<=THRESH_METHOD{2});
    else
        numaps = [];
    end
    %--3
    if ~isempty(numaps)
        disp('number of significant voxels:');
        sum( numaps ~= 0 ),
        
        TMPVOL = zeros( [size(out_analysis.submask), size(tmaps,2)] );
        for i=1:size(tmaps,2)
            tmp = out_analysis.submask;
            tmp(tmp>0) = numaps(:,i);
            TMPVOL(:,:,:,i) = tmp;
        end
        nii=MB;
        nii.img = TMPVOL;
        nii.hdr.dime.datatype = 16;
        nii.hdr.hist = MB.hdr.hist;
        nii.hdr.dime.dim(5) = size(tmaps,2);
        save_untouch_niiz(nii,[out_folder_name,'/',vfield,'_',THRESH_METHOD{1},'.nii']); 
    end

elseif strcmpi( param_type,'mat2d' )

    %--1
    for i=1:size(tmaps,2)
        tmap2d = reshape(tmaps(:,i),matdims);
        writematrix(tmap2d,[out_folder_name,'/',vfield,'_2d_',num2str(i),'_unthresh.txt']); 
    end

    %--2
    for i=1:size(tmaps,2)
        pmap2d = reshape(pmaps(:,i),matdims);
        writematrix(pmap2d,[out_folder_name,'/',pfield,'_2d_',num2str(i),'_unthresh.txt']); 
    end
    
    if strcmpi(THRESH_METHOD{1},'FDR')
        [~,th] = fdr( pmaps,'p',THRESH_METHOD{2},0 );
        numaps = tmaps .* th;
    elseif strcmpi(THRESH_METHOD{1},'CLUST')
        error('cannot spatially cluster conn matrices')
    elseif strcmpi(THRESH_METHOD{1},'UNCORR')
        numaps = tmaps .* double(pmaps<=THRESH_METHOD{2});
    else
        numaps = [];
    end
    %--3
    if ~isempty(numaps)
        disp('number of significant voxels:');
        sum( numaps ~= 0 ),
        tmaps_thresh = numaps;
        
        for i=1:size(numaps,2)
            tmap2d_thresh = reshape(tmaps_thresh(:,i),matdims);
            writematrix(tmap2d_thresh,[out_folder_name,'/',vfield,'_2d_thresh_',num2str(i),'_',THRESH_METHOD{1},'.txt']); 
        end
    end
end

%% storing scores...

    score_arr = NaN*ones( numel(subject_list), 2*size(tmaps,2) );

    if strcmpi(THRESH_METHOD{1},'FDR')
        [~,th] = fdr( pmaps,'p',THRESH_METHOD{2},0 );
        numaps = tmaps .* th;
    elseif strcmpi(THRESH_METHOD{1},'CLUST')
        for i=1:size(tmaps,2)
           tmp = out_analysis.submask;
           tmp(tmp>0) = double(pmaps(:,i)<=0.005) .* tmaps(:,i);
           tmp = clust_up( tmp, THRESH_METHOD{2} );
           numaps(:,i) = tmp(out_analysis.submask>0);
        end
    elseif strcmpi(THRESH_METHOD{1},'UNCORR')
        numaps = tmaps .* double(pmaps<=THRESH_METHOD{2});
    else
        numaps = [];
    end
    %--3
    xn=[];
    xp=[];
    if ~isempty(numaps)
        for i=1:size(tmaps,2)
            if sum(numaps(:,i)<0)>1
                score_arr( :, 2*(i-1)+1 ) = mean( datamat(numaps(:,i)<0,:),1 );
                xn= mean( D2(numaps(:,i)<0,:),1 );
            else
                xn=[];
            end
            if sum(numaps(:,i)>0)>1
                score_arr( :, 2*i ) = mean( datamat(numaps(:,i)>0,:),1 );
                xp= mean( D2(numaps(:,i)>0,:),1 );
            else
                xp=[];
            end

            %--- plotting stage: to be augmented later ---%
            figure, 
            if ~isempty(xn)
            subplot(1,2,1); hold on; title(sprintf('x(%u), negative effects',i));
            plot(X2(:,i), xn,'ok');
            end
            if ~isempty(xp)
            subplot(1,2,2); hold on; title(sprintf('x(%u), positive effects',i));
            plot(X2(:,i), xp,'ok');
            end
        end

    end
    writematrix(score_arr,[out_folder_name,'/score_array_thresh_',num2str(i),'_',THRESH_METHOD{1},'.txt']); 
