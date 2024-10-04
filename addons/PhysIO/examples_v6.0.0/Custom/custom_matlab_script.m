%% Example script using PhysIO with Matlab only (no SPM needed)
%  For documentation of the parameters, see also tapas_physio_new (e.g., via edit tapas_physio_new)
%path: /Users/tomschweizer/Documents/MATLAB/tapas-master/examples_v6.0.0/PhysIO/Custom
%Written by Andrew Hadfield 2024 to evaluate physIO for a task based
%paradigm

%% Create default parameter structure with all fields
physio = tapas_physio_new();


%% Iterate through Individual Particpants
e = dir('/Users/tomschweizer/Documents/MATLAB/tapas-master/examples_v6.0.0/PhysIO/Custom/logfiles/sub*');
puls = dir('/Users/tomschweizer/Documents/MATLAB/tapas-master/examples_v6.0.0/PhysIO/Custom/logfiles/*.puls.log');
resp = dir('/Users/tomschweizer/Documents/MATLAB/tapas-master/examples_v6.0.0/PhysIO/Custom/logfiles/*.resp.log');
folder = '/Users/tomschweizer/Documents/MATLAB/tapas-master/examples_v6.0.0/PhysIO/Custom/logfiles/';

for i = 1:numel(e)
sub = e(i).name(1:8); % participant and session identifiers for respiratory and pulse data
ses = e(i).name(10:18);
if contains(sub, 'sub-092C') %exclusions
else
M=load_untouch_niiz('/Users/tomschweizer/Desktop/TBIfMRI/ATH_BH/fmri_proc/_group_level/masks/pipe_Base1/func_brain_mask_grp.nii');
mask = double(M.img);
load(['/Users/tomschweizer/Desktop/TBIfMRI/ATH_BH/fmri_proc/',sub '_' ses,'/func_proc_p2/pipe_Base1/func1_fullproc.mat']);

    close all;


    combined = e(i).name(1:18);
    initnamepuls = ['/Users/tomschweizer/Documents/MATLAB/tapas-master/examples_v6.0.0/PhysIO/Custom/participants/' sub '_' ses '_bold_BH.puls.1D'];
    initnameresp = ['/Users/tomschweizer/Documents/MATLAB/tapas-master/examples_v6.0.0/PhysIO/Custom/participants/' sub '_' ses '_bold_BH.resp.1D'];
    currentpuls = ['/Users/tomschweizer/Documents/MATLAB/tapas-master/examples_v6.0.0/PhysIO/Custom/logfiles/',sub '_' ses '_bold_BH.puls.log'];
    currentresp = ['/Users/tomschweizer/Documents/MATLAB/tapas-master/examples_v6.0.0/PhysIO/Custom/logfiles/',sub '_' ses '_bold_BH.resp.log'];
    folderPath = fullfile('/Users/tomschweizer/Documents/MATLAB/tapas-master/examples_v6.0.0/PhysIO/Custom/outputs', combined);
    mkdir(folderPath);

    
    physio = tapas_physio_new();
    %Getting errors with this. After break, automate this in order to do
    %batch processing of all log files. Should generate regressors and
    %figures for each one, and then let us check how much global signal is
    %represented by the regressors


%% Individual Parameter settings. Modify to your need and remove default settings
%Anything that is commented out was part of the original code base but
%didn't make the final draft of the pipeline
physio.save_dir = {folderPath}; %'physio_out'
physio.log_files.vendor = 'Custom';
physio.log_files.cardiac = {currentpuls}; %'sub-035C_ses-2-bos_bold_BH.puls.log'
physio.log_files.respiration = {currentresp};
physio.log_files.scan_timing = {''}; % sample rate is 400 but this is expecting log file ,6*60*400, samples for each 
physio.log_files.sampling_interval = 0.0025;
physio.log_files.relative_start_acquisition = 0;
physio.log_files.align_scan = 'last';
physio.scan_timing.sqpar.Nslices = 32; %32, number of axial slices, inside fmRI files 
physio.scan_timing.sqpar.TR = 2;
physio.scan_timing.sqpar.Ndummies = 0;
physio.scan_timing.sqpar.Nscans = size(volmatF,2)+3; %check within fMRI in MRICRON, could be different, may change per person
physio.scan_timing.sqpar.onset_slice = 1;
physio.scan_timing.sync.method = 'nominal';
%physio.scan_timing.sync.method = 'scan_timing_log'; %%not sure on this one, check this
physio.preproc.cardiac.modality = 'PPU';
physio.preproc.cardiac.filter.include = false;
physio.preproc.cardiac.filter.type = 'butter';
physio.preproc.cardiac.filter.passband = [0.3 9];
physio.preproc.cardiac.initial_cpulse_select.method = 'load_from_logfile'; 
% physio.preproc.cardiac.initial_cpulse_select.max_heart_rate_bpm = 90; %might be prudent to change this for breath hold stimuli?
% physio.preproc.cardiac.initial_cpulse_select.file = 'initial_cpulse_kRpeakfile.mat'; %only used if method is manual or load
% physio.preproc.cardiac.initial_cpulse_select.min = 0.4; %could need to adjust for breath hold
% physio.preproc.cardiac.posthoc_cpulse_select.method = 'manual'; %off if manual peak selection not needed
% physio.preproc.cardiac.posthoc_cpulse_select.file = 'manualpeaks.mat'; %added this in
% physio.preproc.cardiac.posthoc_cpulse_select.percentile = 80; %possibly need to adjust these 3 to account for breath hold
% physio.preproc.cardiac.posthoc_cpulse_select.upper_thresh = 60;
% physio.preproc.cardiac.posthoc_cpulse_select.lower_thresh = 60;
physio.preproc.respiratory.filter.passband = [0.01 2];
physio.preproc.respiratory.despike = true; %do we want remove spikes using a sliding window?
physio.model.orthogonalise = 'none'; %can include specific physio noise models
% physio.model.censor_unreliable_recording_intervals = false; %unsure here
physio.model.output_multiple_regressors = 'multiple_regressors.txt';
physio.model.output_physio = 'physio.mat';
physio.model.retroicor.include = true;
physio.model.retroicor.order.c = 3;
physio.model.retroicor.order.r = 4;
physio.model.retroicor.order.cr = 1;
physio.model.rvt.include = true; %should we include this off the start?
physio.model.rvt.method = 'hilbert';
physio.model.rvt.delays = 0;
physio.model.hrv.include = true; %maybe include this after a run of RVT?
physio.model.hrv.delays = 0;
% physio.model.noise_rois.include = false; %havent used this model before?
% physio.model.noise_rois.thresholds = 0.9;
% physio.model.noise_rois.n_voxel_crop = 0;
% physio.model.noise_rois.n_components = 1;
% physio.model.noise_rois.force_coregister = 1;
% physio.model.movement.include = false; %if we work with unprocessed data, may need to include this
% physio.model.movement.order = 6;
% physio.model.movement.censoring_threshold = 0.5;
% physio.model.movement.censoring_method = 'FD';
% physio.model.other.include = false;
physio.verbose.level = 1; %0 if no plot, 1 if just the most important, 2 if you want a bunch, 3 is all
physio.verbose.process_log = cell(0, 1);
physio.verbose.fig_handles = zeros(1, 0);
physio.verbose.use_tabs = false;
physio.verbose.show_figs = false; %was false
physio.verbose.save_figs = false;
physio.verbose.close_figs = false;
physio.ons_secs.c_scaling = 1;
physio.ons_secs.r_scaling = 1;
physio.version = 'R2023a-v9.14.0.2337262'; %changed this to my current matlab version. Was 2022a-v8.1.0

%% Run physiological recording preprocessing and noise modeling
physio = tapas_physio_main_create_regressors(physio);
end
end



% 
% xdes = zeros(30,2);
% xdes(  1: 18,1) = 1;
% xdes( 19: 26,2) = 1;
% xdes = repmat(xdes,6,1);
% xhrf = design_to_hrf( xdes, 2.0, [5 15]);
% xhrf = [xhrf;  zeros(size(volmatF,2) - size(xdes,1) +3,2)  ];
% 
% Xsignal = [ xhrf(4:end,:), physio.model.R(4:end,20) ];
% 
% 
% 
% 
% out = GLM_model_fmri( volmatF, 0, [], Xsignal );
% 
% axial_plot( out.Tmap_signl(:,1), mask, 6, 1, 1); colormap jet; colorbar;
% axial_plot( out.Tmap_signl(:,2), mask, 6, 1, 1); colormap jet; colorbar;
% axial_plot( out.Tmap_signl(:,3), mask, 6, 1, 1); colormap jet; colorbar;
% 
% return;
% 
% [u,l,v] = svd( zscore( volmatF')' ,'econ' );
% figure,plot( zscore([physio.model.R(4:end,20) v(:,1)]),'.-' );
% legend('rrf','bold');
% corr(physio.model.R(4:end,20), v(:,1)*sign(mean(u(:,1)))),

