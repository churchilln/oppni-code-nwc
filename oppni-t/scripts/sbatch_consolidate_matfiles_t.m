close all;
clear;
disp('go onn');

e = dir('texturial_data/sub-*t1.nii');

% divvying up the roi indices
fullixes = [10:13 17:18, 49:54, 11101:11141 11143:11175, 12101:12141 12143:12175];
VBMAX = ceil( numel(fullixes)/5 );

fid = fopen('sbatch_consolidate_report_t.txt','w');

for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        e(i).name,
        % subject-specific information...
        prfix = e(i).name(1:21);    
        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_optim_t/%s',prfix);

        f = dir(sprintf('%s_vbatch*.mat',outfile));

        if isempty(f)
            f2 = dir(sprintf('%s*.mat',outfile));
            if isempty(f2)
                fprintf(fid,'%u. %s : (0) files found\n',i,prfix);
            else
                fprintf(fid,'%u. %s : (1) complete file found\n',i,prfix);
            end
        else
            n = numel(f);
            if n==VBMAX
                % consolidate
                fprintf(fid,'%u. %s : (all) file segments found (of %u) consolidating...!\n',i,prfix,VBMAX);
                %
                outtmp.metrics_av = [];
                outtmp.metrics_sd = [];
                outtmp.opt_stats = [];
                outtmp.CVfull = [];
                outtmp.SZfull = [];
                outtmp.sizecheck = [];
                %
                for vb = 1:VBMAX
                    if ~exist(sprintf('%s_vbatch%uof%u_Output_Roi_byparc.mat',outfile,vb,VBMAX),'file')
                        error('seomthing went wrong! expected to find "%s_vbatch%uof%u_Output_Roi_byparc.mat" but didnt.\n',outfile,vb,VBMAX)
                    end
                    load(sprintf('%s_vbatch%uof%u_Output_Roi_byparc.mat',outfile,vb,VBMAX));
                    outtmp.metrics_av = [outtmp.metrics_av; out.metrics_av];
                    outtmp.metrics_sd = [outtmp.metrics_sd; out.metrics_sd];
                    outtmp.opt_stats  = [outtmp.opt_stats; out.opt_stats];
                    outtmp.CVfull = cat(3,outtmp.CVfull,out.CVfull);
                    outtmp.SZfull = cat(3,outtmp.SZfull,out.SZfull);
                    outtmp.sizecheck = [outtmp.sizecheck; out.sizecheck];
                end
                out = outtmp; clear outtmp;
                save(sprintf('%s_Output_Roi_byparc.mat.mat',outfile),'out');
                if exist(sprintf('%s_Output_Roi_byparc.mat.mat',outfile),'file')
                    for vb = 1:VBMAX
                        unix(sprintf('rm %s_vbatch%uof%u_Output_Roi_byparc.mat',outfile,vb,VBMAX));
                    end
                end
            else
                fprintf(fid,'%u. %s : (%u) file segments found (of %u)\n',i,prfix,n,VBMAX);
            end
        end

    else
        disp('not on list!');
    end
end
fclose(fid);
disp('donne');
