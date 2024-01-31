function P4_fmri_group_level_analysis( inputfile, pipelinefile, paramlist, outpath, param, volno, design_mat, analysis_model, model_contrast, THRESH_METHOD, out_folder_name, censor )
%
% Input:
%      inputfile : string, giving name of input textfile listing data to process
%      pipelinefiles : string, giving name of pipeline textfile specifying which pipeline steps to apply to the data
%      paramlist : string, giving name of parameter file specifying special arguments to use in pipelines
%      outpath : string, specifying destination directory for processed outputs
%
%            volno: index number among contrasts to analyze, if multiple ones
%                           or contrast vector
%       design_mat: design matrix in table format; should have header w variable names
%       analysis_model: string specifying analysis model ('Ttest' or 'GLM')
%      model_contrast: string specifying contrast, based on design matrix fields
%                   (default -- none, do a 1-sample contrast)
%       THRESH_METHOD: threshold method -- 2-element cell array specifying type, critical value}
%                          can be uncorrected: {'UNCORR',[p-value]}
%                               fdr-corrected: {'FDR',[q-value]}
%                   or cluster-size corrected: {'CLUST',[minimum cluster size]}
%                   (default -- none)
%         out_folder_name: string specifying prefix for output files
%                  format is ['Analysis_',<outname>]


% initializing structure and running checkes
[subject_list, ~, PipeStruct_aug, ParamStruct_aug] = P0_fmri_populateDirectories( inputfile, pipelinefile, paramlist, outpath );

if isempty(censor)
    censor = zeros(numel(subject_list),1);
end

% now augmenting outpath... do this after P0!
outpath = fullfile(outpath,'fmri_proc'); % subdir should be fmri_proc

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
            volno = regexp( volno,'-','split');
            bb = str2num( strrep(volno{1},'+',',') );
            aa = str2num( strrep(volno{2},'+',',') );
            datamat(:,i) = mean(x.out_analysis.image.(param)(:,bb),2) - mean(x.out_analysis.image.(param)(:,aa),2);
        end
    elseif strcmpi( param_type,'mat2d' )
        datamat(:,i) = reshape( x.out_analysis.mat2d.(param), [],1);
    else
        error('??')
    end
end


if isempty(design_mat)
    Xdes = [];
    heads = [];
else
    if ischar(design_mat)
        if exist(design_mat,'file')
            design_mat = readtable(design_mat);
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
end

if strcmpi( analysis_model, 'GLM' )
    
    if isempty(model_contrast) error('need to specify a model contrast for GLMs!'); end
    if isempty(Xdes) error('need to load a design matrix for GLMs!'); end
    
    for i=1:numel(model_contrast)
        ix = find( strcmpi( model_contrast{i}, heads) );
        if isempty(ix) error('cannot find variable "%s" in design matrix!\n',model_contrast{i}); end
        Xdes_new(:,i) = Xdes(:,ix);
    end

    D = datamat(:,censor==0);
    % check for bad/missing values
    ixdrop_D = mean(~isfinite(D),1)>0.10;
    fprintf('number of volumes with more than 10% missing data: %s\n',numel(ixdrop_D));
    X = Xdes_new(censor==0,:);
    ixdrop_X = ~isfinite(mean(X,2));
    fprintf('number of design rows with missing data: %s\n',numel(ixdrop_X));
    %
    ixdrop_XD = unique( [ixdrop_D(:); ixdrop_X(:)]);
    fprintf('discarding total: %s\n',numel(ixdrop_XD));
    %
    D(:,ixdrop_XD) = [];
    X(ixdrop_XD,:) = [];

    % checking design matrix
    fprintf('matrix condition (unnormed): %f\n',cond(X));
    fprintf('matrix condition (normed): %f\n',cond(X./sqrt(sum(X.^2,1))));
    % checking for outliers, design matrix
    fprintf('skewness of cols:\n');
    skewness(X),
    fprintf('kurtosis of cols:\n')
    kurtosis(X),


    out_analysis = GLM_gl( datamat(:,censor==0), Xdes_new(censor==0,:) );
    
elseif strcmpi( analysis_model, 'Ttest' )
    
    if isempty(model_contrast)
        disp('no contrast. defaulting to 1-sample t-test');
        Xdes_new = 0;
    else
        if numel(model_contrast>0) error('can only have a single (categorical) predictor for T-testing'); end
        ix = find( strcmpi( model_contrast{1}, heads) );
        if isempty(ix) error('cannot find variable "%s" in design matrix!\n',model_contrast{i}); end
        Xdes_new(:,1) = Xdes(:,ix);
        if numel(unique(Xdes_new)) ~= 2
            error('for t-test, predictor must be categorical');
        end
    end

    if numel(Xdes_new)>1
    out_analysis = Ttest_gl( datamat(:,censor==0), Xdes_new(censor==0,:) );
    else
    out_analysis = Ttest_gl( datamat(:,censor==0), Xdes_new );
    end
end

mkdir_r(out_folder_name);
out_analysis.submask = x.out_analysis.submask;
save([out_folder_name,'/out_analysis.mat'],'out_analysis');

%%

if strcmpi( param_type,'image' )
    
    %--1
    TMPVOL = zeros( [size(out_analysis.submask), size(out_analysis.tstat,2)] );
    for i=1:size(out_analysis.tstat,2)
        tmp = out_analysis.submask;
        tmp(tmp>0) = out_analysis.tstat(:,i);
        TMPVOL(:,:,:,i) = tmp;
    end
    nii=MB;
    nii.img = TMPVOL;
    nii.hdr.dime.datatype = 16;
    nii.hdr.hist = MB.hdr.hist;
    nii.hdr.dime.dim(5) = size(out_analysis.tstat,2);
    save_untouch_niiz(nii,[out_folder_name,'/tmaps_unthresh.nii']); 
    
    %--2
    TMPVOL = zeros( [size(out_analysis.submask), size(out_analysis.tstat,2)] );
    for i=1:size(out_analysis.tstat,2)
        tmp = out_analysis.submask;
        tmp(tmp>0) = out_analysis.tstat_p(:,i);
        TMPVOL(:,:,:,i) = tmp;
    end
    nii=MB;
    nii.img = TMPVOL;
    nii.hdr.dime.datatype = 16;
    nii.hdr.hist = MB.hdr.hist;
    nii.hdr.dime.dim(5) = size(out_analysis.tstat,2);
    save_untouch_niiz(nii,[out_folder_name,'/pmaps_unthresh.nii']); 
    
    if strcmpi(THRESH_METHOD{1},'FDR')
        [~,th] = fdr( out_analysis.tstat_p,'p',THRESH_METHOD{2},0 );
        numaps = out_analysis.tstat .* th;
    elseif strcmpi(THRESH_METHOD{1},'CLUST')
        for i=1:size(out_analysis.tstat,2)
           tmp = out_analysis.submask;
           tmp(tmp>0) = double(out_analysis.tstat_p(:,i)<=0.005) .* out_analysis.tstat(:,i);
           tmp = clust_up( tmp, THRESH_METHOD{2} );
           numaps(:,i) = tmp(out_analysis.submask>0);
        end
    elseif strcmpi(THRESH_METHOD{1},'UNCORR')
        numaps = out_analysis.tstat .* double(out_analysis.tstat_p<=THRESH_METHOD{2});
    else
        numaps = [];
    end
    %--3
    if ~isempty(numaps)
        disp('number of significant voxels:');
        sum( numaps ~= 0 ),
        
        TMPVOL = zeros( [size(out_analysis.submask), size(out_analysis.tstat,2)] );
        for i=1:size(out_analysis.tstat,2)
            tmp = out_analysis.submask;
            tmp(tmp>0) = numaps(:,i);
            TMPVOL(:,:,:,i) = tmp;
        end
        nii=MB;
        nii.img = TMPVOL;
        nii.hdr.dime.datatype = 16;
        nii.hdr.hist = MB.hdr.hist;
        nii.hdr.dime.dim(5) = size(out_analysis.tstat,2);
        save_untouch_niiz(nii,[out_folder_name,'/tmaps_',THRESH_METHOD{1},'.nii']); 
    end

elseif strcmpi( param_type,'mat2d' )

    %--1
    tmaps = out_analysis.tstat;
    for i=1:size(tmaps,2)
        tmap2d = reshape(tmaps(:,i),matdims);
        writematrix(tmap2d,[out_folder_name,'/tmap2d_',num2str(i),'_unthresh.txt']); 
    end

    %--2
    pmaps = out_analysis.tstat_p;
    for i=1:size(tmaps,2)
        pmap2d = reshape(pmaps(:,i),matdims);
        writematrix(pmap2d,[out_folder_name,'/pmap2d_',num2str(i),'_unthresh.txt']); 
    end
    
    if strcmpi(THRESH_METHOD{1},'FDR')
        [~,th] = fdr( out_analysis.tstat_p,'p',THRESH_METHOD{2},0 );
        numaps = out_analysis.tstat .* th;
    elseif strcmpi(THRESH_METHOD{1},'CLUST')
        error('cannot spatially cluster conn matrices')
    elseif strcmpi(THRESH_METHOD{1},'UNCORR')
        numaps = out_analysis.tstat .* double(out_analysis.tstat_p<=THRESH_METHOD{2});
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
            writematrix(tmap2d_thresh,[out_folder_name,'/tmap2d_thresh_',num2str(i),'_',THRESH_METHOD{1},'.txt']); 
        end
    end
end

%% storing scores...

    score_arr = NaN*ones( numel(subject_list), 2*size(tmaps,2) );

    if strcmpi(THRESH_METHOD{1},'FDR')
        [~,th] = fdr( out_analysis.tstat_p,'p',THRESH_METHOD{2},0 );
        numaps = out_analysis.tstat .* th;
    elseif strcmpi(THRESH_METHOD{1},'CLUST')
        for i=1:size(out_analysis.tstat,2)
           tmp = out_analysis.submask;
           tmp(tmp>0) = double(out_analysis.tstat_p(:,i)<=0.005) .* out_analysis.tstat(:,i);
           tmp = clust_up( tmp, THRESH_METHOD{2} );
           numaps(:,i) = tmp(out_analysis.submask>0);
        end
    elseif strcmpi(THRESH_METHOD{1},'UNCORR')
        numaps = out_analysis.tstat .* double(out_analysis.tstat_p<=THRESH_METHOD{2});
    else
        numaps = [];
    end
    %--3
    if ~isempty(numaps)
        for i=1:size(numaps,2)
            if sum(numaps(:,i)<0)>1
                score_arr( :, 2*(i-1)+1 ) = mean( datamat(numaps(:,i)<0,:),1 );
            end
            if sum(numaps(:,i)>0)>1
                score_arr( :, 2*i ) = mean( datamat(numaps(:,i)>0,:),1 );
            end
        end
    end
    writematrix(score_arr,[out_folder_name,'/score_array_thresh_',num2str(i),'_',THRESH_METHOD{1},'.txt']); 
