function [Xreg,stat] = rvhr_PI1( Funcfile, odir, pulsfile, respfile, acqpar, ParamCell)
%
% .ricor_PI1:
% .regression of cardiac and respiratory signals' phase effects with PhysIO utilities 
% .uses "tapas_physio_main_create_regressors" standalone function, then applies slicewise regression

if numel(ParamCell)==1 && isempty(ParamCell{1})
    RV_incl = true;
    HR_incl = true;
elseif numel(ParamCell)==2
    
    ixr=find(contains(ParamCell,'RV'));
    ixh=find(contains(ParamCell,'HR'));
    
    if isempty(ixr)
        RV_incl = false;
    else
        RV_incl = true;
    end
    if isempty(ixh)
        HR_incl = false;
    else
        HR_incl = true;
    end
    
else
    error('need to specify 2 fields for RICOR: e.g., CARD-2,RESP-2,INTER-0')
end

pref = [odir,'/__opptmp_p2func_rvhr'];
% build directory struct recursively
unix(sprintf('mkdir -p %s',pref));

if isempty(acqpar.physamp_msec)
    error('This module needs the physio sampling rate (PHYSAMP_MSEC in the input file)!');
end

%---

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
        physio.model.retroicor.include = false;
        physio.model.rvt.include = RV_incl; %should we include this off the start?
        physio.model.rvt.method = 'hilbert';
        physio.model.rvt.delays = 0;
        physio.model.hrv.include = HR_incl; %maybe include this after a run of RVT?
        physio.model.hrv.delays = 0;
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
        % trimmed array of regressors for mid-slice
        x = load(sprintf('%s/physio_slice0%02u.mat',pref, round(Nslc_eff/2) ));
        Xreg = x.physio.model.R( (1+acqpar.ndrop(1)) : (end-acqpar.ndrop(2)) ,:);


    unix(sprintf('rm -rf %s',pref));

%----

% in cases of nullus datum
ixkep = sum(isfinite(Xreg),1)>0.5*size(Xreg,1);
Xreg = Xreg(:,ixkep);

rr = rank(Xreg(:,:));
[u,l,~]=svd(Xreg);

stat = [rr];
