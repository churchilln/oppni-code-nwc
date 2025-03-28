function ricor_PI1( Funcfile, prefix, odir, pulsfile, respfile, acqpar, ParamCell)
%
% .ricor_PI1:
% .regression of cardiac and respiratory signals' phase effects with PhysIO utilities 
% .uses "tapas_physio_main_create_regressors" standalone function, then applies slicewise regression

if numel(ParamCell)==1 && isempty(ParamCell{1})
    CARD_ord  = 2;
    RESP_ord  = 2;
    INTER_ord = 0;
elseif numel(ParamCell)==3
    
    ixc=find(contains(ParamCell,'CARD-'));
    ixr=find(contains(ParamCell,'RESP-'));
    ixi=find(contains(ParamCell,'INTER-'));
    
    if isempty(ixc) || isempty(ixr) || isempty(ixi)
        error('something missing in RICOR params');
    else
        CARD_ord = str2num(ParamCell{ixc}(6:end));
        RESP_ord = str2num(ParamCell{ixr}(6:end));
        INTER_ord = str2num(ParamCell{ixi}(6:end));
    end    
else
    error('need to specify 3 fields for RICOR: e.g., CARD-2,RESP-2,INTER-0')
end

pref = [odir,'/__opptmp_p2func_ricor'];

spec_case = {'alt+z','alt+z2','alt-z','alt-z2','seq+z','seq-z'};

if isempty(acqpar.physamp_msec)
    error('This module needs the physio sampling rate (PHYSAMP_MSEC in the input file)!');
end

if ~exist(sprintf('%s/%s_ricor.nii.gz',odir,prefix),'file')

    % build directory struct recursively
    unix(sprintf('mkdir -p %s',pref));

    % copyover func file
        unix(sprintf('cp %s %s/func_unric.nii.gz',Funcfile,pref))
        V=load_untouch_niiz(sprintf('%s/func_unric.nii.gz',pref));
        fvol = double(V.img);
        Nslc = size(fvol,3);
        Nt   = size(fvol,4);

    % get "effective" number of slices and their timing -- allows for multiband
        minTE = acqpar.tr_msec/Nslc; % smallest time interval between TEs, if not multiband
        if sum(strcmpi(acqpar.tpatt,spec_case))>0
            % effective "number of slices" for interpolation - same as actual number of slices 
            Nslc_eff = Nslc;
        else
            clear TEsum NVl;
            a = sort(acqpar.tpatt); % sorted list of delays
            ig=1;
            TEsum(ig) = a(1); % sum of delay values @ this interval
            NVl(ig) = 1;     %  number of delay values @ this interval
            for i=2:numel(a)
                % start a new interval if gap is big enough from previous
                if abs(a(i)-a(i-1)) > 0.95*minTE
                    ig=ig+1;
                    TEsum(ig) = a(i);
                    NVl(ig) = 1;
                % otherwise increment for averaging later
                else
                    TEsum(ig)=TEsum(ig)+a(i);
                    NVl(ig)=NVl(ig)+1;
                end
            end
            % get the average
            TEav=TEsum./NVl;
            % effective "number of slices" for interpolation
            Nslc_eff = numel(TEav); 
        end

    % run physio-proc & export regressor files
        puls = load(pulsfile);
        resp = load(respfile);
        % proc -turn 5's into 1's to fit TAPAS guidelines /  insert a column of zeros before the original column
        puls = puls/5;
        puls = [zeros(size(puls, 1), 1), puls];
        % save modified matrix as a .log file
        currentpuls =  sprintf('%s/func.puls.log',pref);
        currentresp =  sprintf('%s/func.resp.log',pref);
        writematrix(puls, currentpuls,'FileType','text');
        writematrix(resp, currentresp,'FileType','text');

        % construct physio structure
        physio = tapas_physio_new();
    
        % Individual Parameter settings...
        physio.save_dir = pref; 
        physio.log_files.vendor = 'Custom';
        physio.log_files.cardiac = {currentpuls}; %'sub-035C_ses-2-bos_bold_BH.puls.log'
        physio.log_files.respiration = {currentresp};
        physio.log_files.scan_timing = {''};
        physio.log_files.sampling_interval = acqpar.physamp_msec./1000; %**
        physio.log_files.relative_start_acquisition = 0;
        physio.log_files.align_scan = 'first';
        physio.scan_timing.sqpar.Nslices = Nslc_eff; %**
        physio.scan_timing.sqpar.TR = acqpar.tr_msec/1000; %**
        physio.scan_timing.sqpar.Ndummies = 0;
        physio.scan_timing.sqpar.Nscans = Nt + sum(acqpar.ndrop); %**
        physio.scan_timing.sqpar.onset_slice = 1:Nslc_eff; %**
        physio.scan_timing.sync.method = 'nominal';
        physio.preproc.cardiac.modality = 'PPU';
        physio.preproc.cardiac.initial_cpulse_select.method = 'load_from_logfile'; 
        physio.preproc.respiratory.filter.passband = [0.01 2];
        physio.preproc.respiratory.despike = true; %do we want remove spikes using a sliding window?
        physio.model.orthogonalise = 'none'; %can include specific physio noise models
        physio.model.output_multiple_regressors = 'multiple_regressors.txt';
        physio.model.output_physio = 'physio.mat';
        physio.model.retroicor.include = true;
        physio.model.retroicor.order.c = CARD_ord;
        physio.model.retroicor.order.r = RESP_ord;
        physio.model.retroicor.order.cr = INTER_ord;
        physio.model.rvt.include = false; %should we include this off the start?
        physio.model.hrv.include = false; %maybe include this after a run of RVT?
        physio.verbose.level = 0; %0 if no plot, 1 if just the most important, 2 if you want a bunch, 3 is all
        physio.verbose.process_log = cell(0, 1);
        physio.verbose.fig_handles = zeros(1, 0);
        physio.verbose.use_tabs = false;
        physio.verbose.show_figs = false; %was false
        physio.verbose.save_figs = false;
        physio.verbose.close_figs = false;
        physio.ons_secs.c_scaling = 1;
        physio.ons_secs.r_scaling = 1;
        
        % Run physiological recording preprocessing and noise modeling
        physio = tapas_physio_main_create_regressors(physio);

        % trimmed array of slice regressors, lag-ordered
        for i=1:Nslc_eff
            x = load(sprintf('%s/physio_slice0%02u.mat',pref,i));
            xmat(:,:,i) = x.physio.model.R( (1+acqpar.ndrop(1)) : (end-acqpar.ndrop(2)) ,:);
        end

        % --slice indexing by acq order
        if sum(strcmpi(acqpar.tpatt,spec_case))>0
            if     strcmpi(acqpar.tpatt,'alt+z')
                sliceidx = [1:2:Nslc, 2:2:Nslc];
            elseif strcmpi(acqpar.tpatt,'alt+z2')
                sliceidx = [2:2:Nslc, 1:2:Nslc];
            elseif strcmpi(acqpar.tpatt,'alt-z')
                if mod(Nslc,2)==0 % even
                    sliceidx = [fliplr(2:2:Nslc) fliplr(1:2:Nslc)];
                else % odd
                    sliceidx = [fliplr(1:2:Nslc) fliplr(2:2:Nslc)];
                end
            elseif strcmpi(acqpar.tpatt,'alt-z2')
                if mod(Nslc,2)==0 % even
                    sliceidx = [fliplr(1:2:Nslc) fliplr(2:2:Nslc)];
                else % odd
                    sliceidx = [fliplr(2:2:Nslc) fliplr(1:2:Nslc)];
                end
            elseif strcmpi(acqpar.tpatt,'seq+z')
                sliceidx = 1:Nslc;
            elseif strcmpi(acqpar.tpatt,'seq-z')
                sliceidx = fliplr( 1:Nslc );
            end
        else
            % if custom, match each slice to the correct delay
            sliceidx = zeros( numel(acqpar.tpatt), 1 );
            for i=1:numel(acqpar.tpatt)
                [~,sliceidx(i)] = min( abs(acqpar.tpatt(i) - TEav) );
            end
            %%deprecated ... assumed all slices were captured uniquely (no MB) 
            %xtt = load(acqpar.tpatt);
            %sliceidx = sortrows([(1:Nslc)' xtt(:)],2,'ascend');
        end

        % slice-by-slice regression, as specificed by tpattern
        % slice idx(i) matched to ith lag-file
        fvol_den = zeros(size(fvol));
        for i=1:Nslc
            i,
            yslc = reshape( fvol(:,:,sliceidx(i),:), [], Nt); %vox x time
            % regress it
            xreg = [ones(Nt,1), xmat(:,:,i) ];
            beta = yslc * xreg /( xreg'*xreg );
            yden = yslc - (beta(:,2:end) * xreg(:,2:end)');
            fvol_den(:,:,sliceidx(i),:) = reshape( yden, size(fvol(:,:,i,:)));
        end

        V.img = fvol_den;
        save_untouch_niiz(V,sprintf('%s/%s_ricor.nii.gz',odir,prefix));

    unix(sprintf('rm -rf %s',pref));
else
    disp('afni-ricor already exists!')
end