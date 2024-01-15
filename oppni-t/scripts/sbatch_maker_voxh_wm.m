close all;
clear;
disp('go onn'); % batching for fixed BW values

e = dir('texturial_data/sub-*t1.nii');

fido = fopen('super_slurm_default.sh','w');

for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        filename = sprintf('tex_sbatch_s%u_default.sl',i);
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
        prfix = e(i).name(1:21);    
        t1_file = sprintf('texturial_data/%s_t1.nii',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii',prfix);
        outfile = sprintf('texturial_matfiles_vox/%s_wm',prfix);
        fprintf(fid,'matlab -nodisplay -nojvm -singleCompThread -r "addpath /lustre03/project/6071833/nwc/NIFTI_tools; addpath /lustre03/project/6071833/nwc/tex_code; texmap_wrapper3( ''%s'', ''%s'', [2 41], 1, ''%s'', 7, 7, [], [], ''HDE'', [0.10] ); exit;"\n',...
            t1_file,mskfile,outfile);
        fclose(fid);

        fprintf(fido,'sbatch %s\n',filename);
    else
        disp('not on list!');
    end
end
fclose(fido);
disp('donne');
