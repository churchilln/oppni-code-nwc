function output = ATBX1( input_file, M0_ref, mask_file, lab_method, lab_order, subtract_type, dataInfo, outname )
%
% =========================================================================
% ATBX1:  code to compute pasl perfusion estimates in OPPNI 
% code framework. Adapted from "ASLtbx" package by Ze Wang
% (https://cfn.upenn.edu/~zewang/ASLtbx.php) and the associated paper.
% Any results should be carefully checked against their codeset if possible.
% =========================================================================
%

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name     = 'ATBX1';
output.attributes.model_descript = 'single-compartment model based on ASLTbx package';

if(nargin==0)
    disp('no inputs - returning attributes');
    return;
end
%----------------------- Default Parameter Checks ------------------------%
%-------------------------------------------------------------------------%

V = load_untouch_niiz(input_file);

if(isempty(mask_file))
    % no mask provided
    disp('no mask provided, estimating (a liberal) one from functional data.') 
    disp('Check outputs carefully to make sure it is correct!');
    mask = double(  mean(V.img,4) > ( 0.2*max(reshape(mean(V.img,4),[],1)) )  );
    M    = V;
    M.img= mask;
else
    % load mask
    M = load_untouch_niiz(mask_file);
    mask = double(M.img>0);
end
% check dimensions of mask/vol
voldims = size(V.img);
for(i=1:3) 
    if( voldims(i) ~= size(M.img,i) ) error(['mask and input volume do not match on dim.',num2str(i)]); end 
end
% to 2D matrix
volmat = nifti_to_mat(V,M);

if( ischar(M0_ref) )
    disp('separate M0 file specified');
    M0vol = load_untouch_niiz(M0_ref);
    m0vct = nifti_to_mat(M0vol,M);
elseif( isnumeric(M0_ref) && isfinite(M0_ref) )
    if( M0_ref > size(voldims(4)) ) error('M0_ref index exceeds size of input file'); end
    if( M0_ref < 1 )                error('M0_ref index must be positive integer');   end
    m0vct = volmat(:,M0_ref); %% assign as M0 volume
    volmat(:,M0_ref)=[]; %% drop from run
else
    error('M0_ref is something unexpected?');
end
% average m0 map
m0vct = mean(m0vct,2);

Ntime = size( volmat,2 );
if( rem(Ntime,2) > 0 ) 
    error('tag and controls are not matched? Should be an even number');
end

% auto-checking on label order
delM_init = mean(volmat(:,1:2:end),2) - mean(volmat(:,2:2:end),2);
% based on median perfusion (must be >0)
if( median(delM_init)>0 ) 
    disp('autocheck: label precedes control');
    if ~isempty(lab_order) && ~strcmpi(lab_order,'TC')
        error('signal autocheck doesnt agree with "TC" label (looks like CT)! you should figure out which is right. Halting!')
    end
    lab_ctl =  1; %%label precedes control
    labix   =  1:2:Ntime;
    conix   =  2:2:Ntime;
else
    disp('autocheck: control precedes label');
    if ~isempty(lab_order) && ~strcmpi(lab_order,'CT')
        error('signal autocheck doesnt agree with "CT" label (looks like TC)! you should figure out which is right. Halting!')
    end
    lab_ctl = -1; %%control precedes label
    conix   =  1:2:Ntime;
    labix   =  2:2:Ntime;
end

% currently only admits 2D PASL option
if contains(lab_method,'PASL')
    fprintf('method %s identified ... carry on!\n',lab_method);
else
    error('sorry! method %s not currently recognized for this package ... check back later!\n',lab_method);
end

% Now pre-initializing parameter values for perfusion estimation (all in msec)
%
% mandatory fields -- users must specify, depends on acquisition
par.TR_MSEC       = NaN;% 2500;  % - repetition time. This is default...may need to adjust
par.TE_MSEC       = NaN;%   12;  % - echo time. This is default...may need to adjust
par.TI1_MSEC      = NaN;%  700;  % - Tagging bolus duration in ms, for the QUIPSS II. This is the default value...may need to adjust
par.TI2_MSEC      = NaN;% 1800;  % - second inversion time; delay time for labeled spin to enter the imaging slice. This is the default value...may need to adjust
par.TSL_MSEC      = NaN; % slice spacing
%
% optional fields -- these go into kinetic modelling equations
par.KM_BLOODT1    = 1600;  % - T1 for blood (msec); ref Lu 04 and Cavusoglu 09 MRI. This is B0 dependent! (BloodT1=1200 at 1.5T...originally BloodT1=1490 at 3.0T / Wang 03)
par.KM_BLOODT2S   =   34;  % - T2s for capillary arterial blood
par.KM_KM_PARCOEF = 0.98;  % - blood/tissue water partition coefficient %*100*60;   %0.9 mL/g
par.KM_PDRAT      = 0.98;  % - proton density ratio blood/gm
par.KM_LABEFF     = 0.90;  % - labeling efficiency, 0.9 for PASL
par.KM_QTI        = 0.85;  % - close to unit, and is set to 0.85 in Warmuth 05
%
par_list = fields(par);

% Updating the par list, based on user input
if( (nargin >= 5) && ~isempty(dataInfo) )
    
   % copying over pre-specified values
   for(i=1:length(par_list))
        if(isfield(dataInfo,par_list{i}))
            disp(['updating ',par_list{i}])
            % updating par. structure
            par.(par_list{i}) = dataInfo.(par_list{i});
        end
        if ~isfinite(par.(par_list{i})) % throw error if mandatory field unspecified
            error('Mandatory par field %s not filled in!',par_list{i});
        end
   end
end

% Slice timing array
if  contains(lab_method,'2D')
    if par.TSL_MSEC ==0
        error('2D ASL requires slice times');
    elseif par.TSL_MSEC ==-1 % auto-esitmation
        % offset in acq. per slice (min-TR - labeltime - delaytime)/#slices
        Slicetime   = (par.TR_MSEC - par.TI2_MSEC)/voldims(3);
    else % input value estimation
        Slicetime   = par.TSL_MSEC;
    end
    % get per-slice adjustments
    tmp         = bsxfun(@times, mask, permute(1:voldims(3),[3 1 2]) ); %% scale voxel values by slice order (start at zero)
    slc_timevct = tmp(mask>0) .* Slicetime; %% vectorized, scaled by slice timing
elseif contains(lab_method,'3D')
    if par.TSL_MSEC <=0 % no relative delays
        slc_timevct = 0;        
    else
        error('3D ASL doesnt have slice offsets');
    end
end

%% COMPUTING PERFUSION WEIGHTED IMAGES + BOLD series

if    ( (isnumeric(subtract_type) && (subtract_type==1)) || strcmpi( subtract_type,'simple') )

    % subtract nearest control from each label
    PERF = volmat(:,labix) - volmat(:,conix);
    
elseif( (isnumeric(subtract_type) && (subtract_type==2)) || strcmpi( subtract_type,'surround') )
    
    if    ( lab_ctl>0 ) 
        % avg of matched(subseq) and prev
        con_mat = ( volmat(:,conix) + [volmat(:,conix(1)) volmat(:,conix(1:end-1))] )./2;
        %
    elseif( lab_ctl<0 )
        % avg of matched(prior) and subseq
        con_mat = ( volmat(:,conix) + [volmat(:,conix(2:end)) volmat(:,conix(end))] )./2;
        %
    end
    % subtract from tag
    PERF = volmat(:,labix) - con_mat;
    
elseif( (isnumeric(subtract_type) && (subtract_type==3)) || strcmpi( subtract_type,'sinc') )
    
    % lab_ctl -1  same as  FirstImageType==1 (control before label)
    Timeshift = 0.5;
    for(n=1:round(Ntime/2)) %% step through volumes
       %
       % 6 point sinc interpolation
       if     lab_ctl>0
           idx=n+[-3 -2 -1 0 1 2];
           normloc=3-Timeshift;
       elseif lab_ctl<0
           idx=n+[-2 -1 0 1 2 3];
           normloc=2+Timeshift;
       end
       idx(idx<1)=1;
       idx(idx>round(Ntime/2))=round(Ntime/2);
       con_img=sinc_interpVec( volmat(:,conix(idx)) ,normloc);
       % now subtract interpolated value from tag
       PERF(:,n) = volmat(:,labix(n)) - con_img;
    end
    
else
    error('subtract_type needs to be (1,2,3) or (simple, surround, sinc)');
end

clear delM_init; %% clear initialized

%% NOW QUANTIFICATION OF CBF

% adjusted TI2, based on TI2 (delay time) & slice-specific delay
TI_adj= (par.TI2_MSEC) + slc_timevct;
% equilibrium magnetization of arterial blood
M0a   = (par.KM_PDRAT * exp(-par.TE_MSEC/par.KM_BLOODT2S)) .* m0vct;
% absolute CBF, using the ASLtbx formulation of TCBF (ml/100g/ms) --> (ml/100g/min):
aCBF = 6000*1000 * bsxfun(@rdivide, par.KM_KM_PARCOEF.*PERF , (2*par.KM_LABEFF*par.TI1_MSEC*par.KM_QTI).*M0a.*exp(-TI_adj./par.KM_BLOODT1) );
% uncalibrated CBF --> if 2D, does a single adjustment for slice-dependent transit times
uCBF = bsxfun(@rdivide, PERF, exp(-TI_adj./par.KM_BLOODT1) );

% Storing outputs
output.aCBF = aCBF;
output.uCBF = uCBF;
% averages too
output.aCBF_mean = mean(aCBF,2);
output.uCBF_mean = mean(uCBF,2);
output.M0        = mean(m0vct,2);

if ~isempty(outname)
    
    % predefine header nii
    nii=V;
    nii.hdr.dime.datatype = 16;
    nii.hdr.hist = V.hdr.hist;
    % initialize 4d tmpvol
    TMPVOL = zeros( [size(mask), round(Ntime/2)] ); 

    % --------------------- saving 4d functional data ---------------------
    nii.hdr.dime.dim(5) = size(TMPVOL,4);

    % --tcbf--
    for(t=1:round(Ntime/2)) 
        tmp=mask;tmp(tmp>0)=output.aCBF(:,t); 
        TMPVOL(:,:,:,t) = tmp; 
    end
    nii.img = TMPVOL;
    save_untouch_niiz(nii,[outname,'_CBF.nii.gz']);
    
    % --perfusion-- 
    for(t=1:round(Ntime/2)) 
        tmp=mask;tmp(tmp>0)=output.uCBF(:,t); 
        TMPVOL(:,:,:,t) = tmp; 
    end
    nii.img = TMPVOL;
    save_untouch_niiz(nii,[outname,'_PWI.nii.gz']);

    % --------------------- saving the average volumes ---------------------
    nii.hdr.dime.dim(5) = 1;
    
    % --mean:tcbf--
    TMPVOL=mask;TMPVOL(TMPVOL>0)=output.aCBF_mean;
    nii.img = TMPVOL;
    save_untouch_niiz(nii,[outname,'_CBF_avg.nii.gz']);
    nii.img = TMPVOL;
    
    % --mean:perfusion--
    TMPVOL=mask;TMPVOL(TMPVOL>0)=output.uCBF_mean;
    nii.img = TMPVOL;
    save_untouch_niiz(nii,[outname,'_PWI_avg.nii.gz']);
end

%%
function y = sinc_interpVec(x,u)
% sinc interpolation function, the input can be a data vector or matrix,
% the number of rows is assumed to be the number of data dimension, the
% number of columns is assumed to be the time points. The interpolation is
% applied column wise.
% Ze Wang @ 6-09-2004 Upenn
[dim,len]=size(x);
[dim2,ulen]=size(u);
if dim2==1
    u=repmat(u,dim,1);
else
    if dim2~=dim disp('We can'' figure out what you want to do\n');return; end;
end
m = 0:len-1;
m=repmat(m,[dim,1]);
for i=1:ulen
    weight=sinc(m- repmat(u(:,i),1,len));
    swei=sum(weight,2);
    if abs(swei(1)-1)>0.1
        weight=weight./repmat(swei,1,len);
    end
  
  y(:,i) = sum(x.*weight, 2);
end
