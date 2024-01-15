close all;
clear;
disp('go onn'); % batching for fixed BW values

mkdir('texturial_matfiles_t_aD');
e = dir('texturial_data/sub-*t1.nii');
fido = fopen('super_slurm_t_ad.sh','w');
load('tex_code/bwarr_t.mat');

kq=0;
vcount=0;
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')

        kq=kq+1;
        prfix = e(i).name(1:21);  
        t1_file = sprintf('texturial_data/%s_t1.nii',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii',prfix);

        outfile = sprintf('texturial_matfiles_t_aD/%s_aD',prfix);
        outfile_fullname = sprintf('%s_Output_Roi_byparc.mat',outfile);
        if ~exist(outfile_fullname,'file')

            vcount = vcount+1;
            filename = sprintf('tex_sbatch_s%u_ad.sl',i);
            fid = fopen(filename,'w');
            % header prep stuff
            fprintf(fid,'#!/bin/bash -l\n');
            fprintf(fid,'#SBATCH --job-name=matlab_test\n');
            fprintf(fid,'#SBATCH --account=def-nwc\n'); % adjust this to match the accounting group you are using to submit jobs
            fprintf(fid,'#SBATCH --time=0-03:00\n');         % adjust this to match the walltime of your job
            fprintf(fid,'#SBATCH --nodes=1\n');      
            fprintf(fid,'#SBATCH --ntasks=1\n');
            fprintf(fid,'#SBATCH --cpus-per-task=1\n');      % adjust this if you are using parallel commands
            fprintf(fid,'#SBATCH --mem=16G\n');             % adjust this according to the memory requirement per node you need
            fprintf(fid,'\n\n');
            % module loadink
            fprintf(fid,'module load StdEnv/2020\n');
            fprintf(fid,'module load matlab/2022a\n');
            fprintf(fid,'\n\n');
            % actual command to execute
            if ~strcmpi( prfix, prfix_arr{kq} )
                error('matchnigueq oafd idsas');            
            end
            fprintf(fid,'matlab -nodisplay -nojvm -singleCompThread -r "addpath /lustre03/project/6071833/nwc/NIFTI_tools; addpath /lustre03/project/6071833/nwc/tex_code; load(''tex_code/bwarr_t.mat''); szoptim = mean(szarr,3); texmap_wrapper3( ''%s'', ''%s'', [10:13 17:18, 49:54,  11101:11141 11143:11175, 12101:12141 12143:12175], 0, ''%s'', [], [], [], [], ''KDE-GT'', [szoptim(:,%u)] ); exit;"\n',...
                t1_file,mskfile,outfile,kq);
            fclose(fid);
    
            fprintf(fido,'sbatch %s\n',filename);
        else
            fprintf('skipping "%s", file already exists\n\n',outfile_fullname);
        end

        outfile = sprintf('texturial_matfiles_t_aD/%s_aDaS',prfix);
        outfile_fullname = sprintf('%s_Output_Roi_byparc.mat',outfile);
        if ~exist(outfile_fullname,'file')

            vcount = vcount+1;
            filename = sprintf('tex_sbatch_s%u_adas.sl',i);
            fid = fopen(filename,'w');
            % header prep stuff
            fprintf(fid,'#!/bin/bash -l\n');
            fprintf(fid,'#SBATCH --job-name=matlab_test\n');
            fprintf(fid,'#SBATCH --account=def-nwc\n'); % adjust this to match the accounting group you are using to submit jobs
            fprintf(fid,'#SBATCH --time=0-03:00\n');         % adjust this to match the walltime of your job
            fprintf(fid,'#SBATCH --nodes=1\n');      
            fprintf(fid,'#SBATCH --ntasks=1\n');
            fprintf(fid,'#SBATCH --cpus-per-task=1\n');      % adjust this if you are using parallel commands
            fprintf(fid,'#SBATCH --mem=16G\n');             % adjust this according to the memory requirement per node you need
            fprintf(fid,'\n\n');
            % module loadink
            fprintf(fid,'module load StdEnv/2020\n');
            fprintf(fid,'module load matlab/2022a\n');
            fprintf(fid,'\n\n');
            % actual command to execute
            if ~strcmpi( prfix, prfix_arr{kq} )
                error('matchnigueq oafd idsas');            
            end
            fprintf(fid,'matlab -nodisplay -nojvm -singleCompThread -r "addpath /lustre03/project/6071833/nwc/NIFTI_tools; addpath /lustre03/project/6071833/nwc/tex_code; load(''tex_code/bwarr_t.mat''); szoptim = median(mean(szarr,3),2); texmap_wrapper3( ''%s'', ''%s'', [10:13 17:18, 49:54,  11101:11141 11143:11175, 12101:12141 12143:12175], 0, ''%s'', [], [], [], [], ''KDE-GT'', [szoptim(:)] ); exit;"\n',...
                t1_file,mskfile,outfile);
            fclose(fid);
    
            fprintf(fido,'sbatch %s\n',filename);
        else
            fprintf('skipping "%s", file already exists\n\n',outfile_fullname);
        end

        outfile = sprintf('texturial_matfiles_t_aD/%s_aDaSaR',prfix);
        outfile_fullname = sprintf('%s_Output_Roi_byparc.mat',outfile);
        if ~exist(outfile_fullname,'file')

            vcount = vcount+1;
            filename = sprintf('tex_sbatch_s%u_adasar.sl',i);
            fid = fopen(filename,'w');
            % header prep stuff
            fprintf(fid,'#!/bin/bash -l\n');
            fprintf(fid,'#SBATCH --job-name=matlab_test\n');
            fprintf(fid,'#SBATCH --account=def-nwc\n'); % adjust this to match the accounting group you are using to submit jobs
            fprintf(fid,'#SBATCH --time=0-03:00\n');         % adjust this to match the walltime of your job
            fprintf(fid,'#SBATCH --nodes=1\n');      
            fprintf(fid,'#SBATCH --ntasks=1\n');
            fprintf(fid,'#SBATCH --cpus-per-task=1\n');      % adjust this if you are using parallel commands
            fprintf(fid,'#SBATCH --mem=16G\n');             % adjust this according to the memory requirement per node you need
            fprintf(fid,'\n\n');
            % module loadink
            fprintf(fid,'module load StdEnv/2020\n');
            fprintf(fid,'module load matlab/2022a\n');
            fprintf(fid,'\n\n');
            % actual command to execute
            if ~strcmpi( prfix, prfix_arr{kq} )
                error('matchnigueq oafd idsas');            
            end
            fprintf(fid,'matlab -nodisplay -nojvm -singleCompThread -r "addpath /lustre03/project/6071833/nwc/NIFTI_tools; addpath /lustre03/project/6071833/nwc/tex_code; load(''tex_code/bwarr_t.mat''); szoptim = median(mean(szarr,3),2); texmap_wrapper3( ''%s'', ''%s'', [10:13 17:18, 49:54,  11101:11141 11143:11175, 12101:12141 12143:12175], 0, ''%s'', [], [], [], [], ''KDE-GT'', [szoptim] ); exit;"\n',...
                t1_file,mskfile,outfile);
            fclose(fid);
    
            fprintf(fido,'sbatch %s\n',filename);
        else
            fprintf('skipping "%s", file already exists\n\n',outfile_fullname);
        end

    else
        disp('not on list!');
    end
end
fclose(fido);
if vcount==0
    disp('no missing files - removing jobfile');
    unix(sprintf('rm super_slurm_ad.sh'))
else
    fprintf('found %u missing files! to run!',vcount);
end
disp('donne');
