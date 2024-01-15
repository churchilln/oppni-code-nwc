function [QCtable_compat QCtable_quant] = QCstats_raw( filepaths, prefixes, outpath, type, lean, minqc )
%
% script collects basic quality control information for minimally
% processed scans, anatomical or functional
% 

[NSmax NRmax] = size(filepaths);

if isempty(prefixes)
    for ns=1:NSmax
        prefixes{ns,1} = sprintf('subj_%u',ns);
    end
end
if nargin<4 || (~strcmpi(type,'anat') && ~strcmpi(type,'func') && ~strcmpi(type,'perf') && ~strcmpi(type,'diff'))
    error('need to specify proper datatype - anat, func, perf or diff!')
end
if nargin<5
    lean=1;
end
if nargin<6
    minqc=0;
end

fprintf('starting %s-raw qc...\n\n',type);

mkdir_r(outpath);
QCtable_compat = table();
if minqc==0
QCtable_quant  = table();
end

kq=0;
for ns=1:NSmax
    
    opath0 = sprintf('%s/%s',outpath,prefixes{ns}); % construct subject directory
    unix(sprintf('mkdir %s',opath0));

    for nr=1:NRmax
    
        if ~isempty(filepaths{ns,nr})

            kq=kq+1; % increment if run exists

            % populating this row
            QCtable_compat.ID{kq}  = prefixes{ns};
            QCtable_quant.ID{kq}   = prefixes{ns};
            QCtable_compat.run(kq) = nr;
            QCtable_quant.run(kq)  = nr;
        
            % port over, minimally proc and check the anatomical scan
            if ~exist(sprintf('%s/img%u.nii',opath0,nr),'file')
                if contains(filepaths{ns},'.nii.gz')
                    unix(sprintf('cp %s %s/img%u.nii.gz',filepaths{ns},opath0,nr))
                    unix(sprintf('gunzip %s/img%u.nii.gz',opath0,nr)); % unzipt
                elseif contains(filepaths{ns},'.nii')
                    unix(sprintf('cp %s %s/img%u.nii',filepaths{ns},opath0,nr))
                else
                    error('unrecognized input file type for:\n\t%s\n',filepaths{ns})
                end
            end

            % ---->>>> basic compatibility checking
        
            % gathering information for qc compatibility stats
            hdr = load_nii_hdr(sprintf('%s/img%u.nii',opath0,nr));
            % now zippit
            unix(sprintf('gzip %s/img%u.nii',opath0,nr));
            % orientation info
            orient_RASO = sign([ hdr.hist.srow_x(1) hdr.hist.srow_y(2) hdr.hist.srow_z(3) ]);
            if    ( prod( orient_RASO == [ 1 1 1] )==1 ) orient_RASO(4) =  1; % RAS=radiologic
            elseif( prod( orient_RASO == [-1 1 1] )==1 ) orient_RASO(4) = -1; % LAS=neurologic
            else orient_RASO(4) = 0; % other/unidentified
            end
            %--- store to qc compatibility struct
            QCtable_compat.xres(kq)  = hdr.dime.pixdim(2);
            QCtable_compat.yres(kq)  = hdr.dime.pixdim(3);
            QCtable_compat.zres(kq)  = hdr.dime.pixdim(4);
            QCtable_compat.tres(kq)  = hdr.dime.pixdim(5);
            %
            QCtable_compat.xdim(kq)  = hdr.dime.dim(2);
            QCtable_compat.ydim(kq)  = hdr.dime.dim(3);
            QCtable_compat.zdim(kq)  = hdr.dime.dim(4);
            QCtable_compat.tdim(kq)  = hdr.dime.dim(5);
            %
            QCtable_compat.bitpx(kq) = hdr.dime.bitpix;
            QCtable_compat.ori_R(kq) = orient_RASO(1);
            QCtable_compat.ori_A(kq) = orient_RASO(2);
            QCtable_compat.ori_S(kq) = orient_RASO(3);
            QCtable_compat.ori_O(kq) = orient_RASO(4);

            % ---->>>> some quantitative checks

            if minqc==0 %%-- if quant checks are requested
        
            if ~exist(sprintf('%s/%s%u_qun.mat',opath0,type,nr),'file')
            
                unix(sprintf('mkdir %s/deriv',opath0));
                if ~exist(sprintf('%s/deriv/img%u_brmask_bin.nii.gz',opath0,nr),'file')
                    % deoblique and skullstrip
                    unix(sprintf('3dWarp -oblique2card -prefix %s/deriv/img%u_deob.nii.gz -cubic %s/img%u.nii.gz', opath0,nr,opath0,nr)); %-wsinc5

                    if strcmpi(type,'anat')
                        unix(sprintf('3dSkullstrip -input %s/deriv/img%u_deob.nii.gz -prefix %s/deriv/img%u_brmask.nii.gz -mask_vol',opath0,nr,opath0,nr));
                        % convert mask to binary, 
                        M0=load_untouch_niiz(sprintf('%s/deriv/img%u_brmask.nii.gz',opath0,nr));
                        M0.img = double(M0.img==6);
                        save_untouch_niiz(M0,sprintf('%s/deriv/img%u_brmask_bin.nii.gz',opath0,nr));
                        unix(sprintf('rm %s/deriv/img%u_brmask.nii.gz',opath0,nr));

                    elseif strcmpi(type,'func')
                        unix(sprintf('3dAutomask -prefix %s/deriv/func%u_brmask_bin.nii.gz %s/deriv/func%u_deob.nii.gz',opath0,nr,opath0,nr));
                        % convert mask to binary, 
                        M0=load_untouch_niiz(sprintf('%s/deriv/func%u_brmask_bin.nii.gz',opath0,nr));
                    end
                else
                    M0=load_untouch_niiz(sprintf('%s/deriv/img%u_brmask_bin.nii.gz',opath0,nr));
                end
                if strcmpi(type,'anat') && ~exist(sprintf('%s/deriv/img%u_seg_GM.nii.gz',opath0,nr),'file') % only available for anat
                
                    V0 = load_untouch_niiz(sprintf('%s/deriv/img%u_deob.nii.gz',opath0,nr));
                    M0 = load_untouch_niiz(sprintf('%s/deriv/img%u_brmask_bin.nii.gz',opath0,nr));
                    V0.img = double(V0.img .* M0.img);
                    save_untouch_niiz(V0,sprintf('%s/deriv/img%u_masked.nii.gz',opath0,nr));
            
                    unix(sprintf('fast -R 0.3 -H 0.1 -t 1 %s/deriv/img%u_masked.nii.gz',opath0,nr));
                
                    % renaming 'em
                    VX = load_untouch_niiz(sprintf('%s/deriv/img%u_masked.nii.gz',opath0,nr));
                    for i=1:3
                        VS = load_untouch_niiz(sprintf('%s/deriv/img%u_masked_pve_%u.nii.gz',opath0,nr,i-1));
                        sval(i,1) = mean(double(VX.img(VS.img>0.5)));
                    end
                    isort = sortrows([(1:3)',sval],2,'ascend');
                    isort = isort(:,1);           % pves indexed by increasing mean intensity
                    tisslist = {'CSF','GM','WM'}; % tissues in increasing order of T1 intensity
                    for i=1:3
                        unix(sprintf('mv %s/deriv/img%u_masked_pve_%u.nii.gz %s/deriv/img%u_seg_%s.nii.gz',opath0,nr, (isort(i)-1), opath0,nr, tisslist{i}));
                    end
                    unix(sprintf('rm %s/deriv/img%u_masked_pve*',opath0,nr));
                end
                %==loadnvolumes==%
                if strcmpi(type,'anat')
                    VG = load_untouch_niiz(sprintf('%s/deriv/img%u_seg_GM.nii.gz',opath0,nr));
                    VW = load_untouch_niiz(sprintf('%s/deriv/img%u_seg_WM.nii.gz',opath0,nr));
                    VC = load_untouch_niiz(sprintf('%s/deriv/img%u_seg_CSF.nii.gz',opath0,nr));
                end
            
                if ~exist(sprintf('%s/deriv/img%u_MASKIN.nii.gz',opath0,nr),'file') || ~exist(sprintf('%s/deriv/img%u_MASKOUT.nii.gz',opath0,nr),'file')
                    % take z-axis cut for signal estimation
                    for z=1:size(M0.img,3) cseca(z,1) = sum(sum( M0.img(:,:,z) )); end
                    [~,izm]=max(cseca);
                    % in-of-brain mask]
                    Mi=M0;
                    Mi.img(:,:,1:izm-1)=0;
                    save_untouch_niiz(Mi,sprintf('%s/deriv/img%u_MASKIN.nii.gz',opath0,nr));
                    % out-of-brain mask
                    ndil = ceil(25/min(M0.hdr.dime.pixdim(2:4)));
                    unix(sprintf('3dmask_tool -dilate_inputs %u -prefix %s/deriv/img%u_brmask_dil.nii.gz -input %s/deriv/img%u_brmask_bin.nii.gz',ndil,opath0,nr,opath0,nr));
                    Mo=load_untouch_niiz(sprintf('%s/deriv/img%u_brmask_dil.nii.gz',opath0,nr));
                    Mo.img = double(Mo.img>0);
                    Mo.img = (1 - Mo.img);
                    Mo.img(:,:,1:izm-1)=0;
                    save_untouch_niiz(Mo,sprintf('%s/deriv/img%u_MASKOUT.nii.gz',opath0,nr));
                    % delete full dilated mask
                    unix(sprintf('rm %s/deriv/img%u_brmask_dil.nii.gz',opath0,nr));
                else
                    Mi=load_untouch_niiz(sprintf('%s/deriv/img%u_MASKIN.nii.gz',opath0,nr));
                    Mo=load_untouch_niiz(sprintf('%s/deriv/img%u_MASKOUT.nii.gz',opath0,nr));
                end
                V0 = load_untouch_niiz(sprintf('%s/deriv/img%u_deob.nii.gz',opath0,nr));
                %
                qun_vol = V0.img;
                qun_msk = M0.img;
                qun_mski= Mi.img;
                qun_msko= Mo.img;
                qun_vdim= V0.hdr.dime.pixdim(2:4);
                if strcmpi(type,'anat')
                    qun_pgm = VG.img;
                    qun_pwm = VW.img;
                    qun_pcf = VC.img;
                    save(sprintf('%s/%s%u_qun.mat',opath0,type,nr),'qun_vol','qun_msk','qun_mski','qun_msko','qun_pgm','qun_pwm','qun_pcf','qun_vdim');
                else
                    save(sprintf('%s/%s%u_qun.mat',opath0,type,nr),'qun_vol','qun_msk','qun_mski','qun_msko','qun_vdim');
                end
            else
                load(sprintf('%s/%s%u_qun.mat',opath0,type,nr));
            end
            if lean>0 % remove derivatives
                unix(sprintf('rm -r %s/deriv',opath0));
            end


            % max-likelihood tissue mask
            if strcmpi(type,'anat')
            [~,tiss_ix] = max( cat(4,qun_pcf,qun_pgm,qun_pwm),[],4);
            end

            % general measures
            for t=1:size(qun_vol,4)
            [QCtable_quant.SNRa(kq)] = mean( double(qun_vol( qun_mski>0 )) ) ./ std( double(qun_vol( qun_msko>0 )) ); % SNR > all-brain mean vs. noise SD
            [QCtable_quant.CoVa(kq)] =  std( double(qun_vol( qun_mski>0 )) ) ./ mean( double(qun_vol( qun_mski>0 )) ); % spatial signal homogeneity
            [QCtable_quant.sACav(kq), QCtable_quant.sACd1(kq), QCtable_quant.sACd2(kq)] = estim_acd_spat( qun_vol, qun_mski ); % signal smoothness
            [QCtable_quant.nACav(kq), QCtable_quant.nACd1(kq), QCtable_quant.nACd2(kq)] = estim_acd_spat( qun_vol, qun_msko ); % noise smoothness
            [QCtable_quant.mvol(kq)] = sum(qun_msk(:)) .* prod(qun_vdim);
            end
           
            if strcmpi(type,'anat')
                % general spatial measures
                [QCtable_quant.SNRa_av(kq)] = mean( double(qun_vol( qun_mski>0 )) ) ./ std( double(qun_vol( qun_msko>0 )) ); % SNR > all-brain mean vs. noise SD
                [QCtable_quant.CoVa_av(kq)] =  std( double(qun_vol( qun_mski>0 )) ) ./ mean( double(qun_vol( qun_mski>0 )) ); % spatial signal homogeneity
                [QCtable_quant.sACav_av(kq), QCtable_quant.sACd1(kq), QCtable_quant.sACd2(kq)] = estim_acd_spat( qun_vol, qun_mski ); % signal smoothness
                [QCtable_quant.nACav_av(kq), QCtable_quant.nACd1(kq), QCtable_quant.nACd2(kq)] = estim_acd_spat( qun_vol, qun_msko ); % noise smoothness
                [QCtable_quant.mvol(kq)] = sum(qun_msk(:)) .* prod(qun_vdim);
                 % anatomical, tissue-specific measures
                [QCtable_quant.SNRc(kq)] = mean( double(qun_vol( tiss_ix==1 )) ) ./ std( double(qun_vol( qun_msko>0 )) ); % SNR > all-brain mean vs. noise SD
                [QCtable_quant.SNRg(kq)] = mean( double(qun_vol( tiss_ix==2 )) ) ./ std( double(qun_vol( qun_msko>0 )) ); % SNR > all-brain mean vs. noise SD
                [QCtable_quant.SNRw(kq)] = mean( double(qun_vol( tiss_ix==3 )) ) ./ std( double(qun_vol( qun_msko>0 )) ); % SNR > all-brain mean vs. noise SD
                [QCtable_quant.CoVc(kq)] =  std( double(qun_vol( tiss_ix==1 )) ) ./ mean( double(qun_vol( tiss_ix==1 )) ); % spatial signal homogeneity
                [QCtable_quant.CoVg(kq)] =  std( double(qun_vol( tiss_ix==2 )) ) ./ mean( double(qun_vol( tiss_ix==2 )) ); % spatial signal homogeneity
                [QCtable_quant.CoVw(kq)] =  std( double(qun_vol( tiss_ix==3 )) ) ./ mean( double(qun_vol( tiss_ix==3 )) ); % spatial signal homogeneity
                [QCtable_quant.CNRcg(kq)] = abs( mean( double(qun_vol( tiss_ix==1 )) ) - mean( double(qun_vol( tiss_ix==2 )) ) ) ./ std( double(qun_vol( qun_msko>0 )) ); % signal quality (MSEratio is very similar / less teails)
                [QCtable_quant.CNRcw(kq)] = abs( mean( double(qun_vol( tiss_ix==1 )) ) - mean( double(qun_vol( tiss_ix==3 )) ) ) ./ std( double(qun_vol( qun_msko>0 )) ); % signal quality (MSEratio is very similar / less teails)
                [QCtable_quant.CNRgw(kq)] = abs( mean( double(qun_vol( tiss_ix==2 )) ) - mean( double(qun_vol( tiss_ix==3 )) ) ) ./ std( double(qun_vol( qun_msko>0 )) ); % signal quality (MSEratio is very similar / less teails)
            elseif strcmpi(type,'func')
                % general spatial measures
                for t=1:size(qun_vol,4)
                    qun_vol_tmp = qun_vol(:,:,:,t);
                    [SNRa(t,1)] = mean( double(qun_vol_tmp( qun_mski>0 )) ) ./ std( double(qun_vol_tmp( qun_msko>0 )) ); % SNR > all-brain mean vs. noise SD
                    [CoVa(t,1)] =  std( double(qun_vol_tmp( qun_mski>0 )) ) ./ mean( double(qun_vol_tmp( qun_mski>0 )) ); % spatial signal homogeneity
                    [sACav(t,1), sACd1(t,1), sACd2(t,1)] = estim_acd_spat( qun_vol_tmp, qun_mski ); % signal smoothness
                    [nACav(t,1), nACd1(t,1), nACd2(t,1)] = estim_acd_spat( qun_vol_tmp, qun_msko ); % noise smoothness
                end
                % time-av
                QCtable_quant.SNRa_tav(kq)  = mean(SNRa);
                QCtable_quant.CoVa_tav(kq)  = mean(CoVa);
                QCtable_quant.sACav_tav(kq) = mean(sACav);
                QCtable_quant.sACd1_tav(kq) = mean(sACd1);
                QCtable_quant.sACd2_tav(kq) = mean(sACd2);
                QCtable_quant.nACav_tav(kq) = mean(nACav);
                QCtable_quant.nACd1_tav(kq) = mean(nACd1);
                QCtable_quant.nACd2_tav(kq) = mean(nACd2);
                % time-vr
                QCtable_quant.SNRa_tsd(kq)  =  std(SNRa);
                QCtable_quant.CoVa_tsd(kq)  =  std(CoVa);
                QCtable_quant.sACav_tsd(kq) =  std(sACav);
                QCtable_quant.sACd1_tsd(kq) =  std(sACd1);
                QCtable_quant.sACd2_tsd(kq) =  std(sACd2);
                QCtable_quant.nACav_tsd(kq) =  std(nACav);
                QCtable_quant.nACd1_tsd(kq) =  std(nACd1);
                QCtable_quant.nACd2_tsd(kq) =  std(nACd2);
                %
                [QCtable_quant.mvol(kq)] = sum(qun_msk(:)) .* prod(qun_vdim);
                % functional, temporal-style measures
                qun_vol_tmp = mean(qun_vol,4)./std(qun_vol,0,4);
                QCtable_quant.SNRt_sav(kq) = mean( qun_vol_tmp(qun_mski>0) );
                QCtable_quant.SNRt_ssd(kq) =  std( qun_vol_tmp(qun_mski>0) );
                %
                [QCtable_quant.tACav(kq), QCtable_quant.tACd1(kq), QCtable_quant.tACd2(kq)] = estim_acd_temp( double(qun_vol_set), double(qun_mski) );
            end
            %-
            clear qun_vol qun_msk qun_mski qun_msko qun_pgm qun_pwm qun_pcf qun_vdim qun_vol_tmp;
            clear SNRa CoVa sACav sACd1 sACd2 nACav nACd1 nACd2;

            end

        end
    end
end

save(sprintf('%s/QCtable_%s_compat.mat',outpath,type),'QCtable_compat')
if minqc==0
save(sprintf('%s/QCtable_%s_quant_raw.mat',outpath,type), 'QCtable_quant' )
end

disp('finished!\n')
