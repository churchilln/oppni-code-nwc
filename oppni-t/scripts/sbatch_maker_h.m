close all;
clear;
disp('go onn'); % batching for optimized BW value - HDE (more time intensive!!)

e = dir('texturial_data/sub-*t1.nii');

% divvying up the roi indices
fullixes = [10:13 17:18, 49:54, 11101:11141 11143:11175, 12101:12141 12143:12175];
VBMAX = ceil( numel(fullixes)/5 );

for i=1:5%numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')

        % subject-specific information...
        prfix = e(i).name(1:21);    
        t1_file = sprintf('texturial_data/%s_t1.nii',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii',prfix);
        outfile = sprintf('texturial_matfiles_optim_h/%s',prfix);

        % master job submitter for this subject
        fido = fopen(sprintf('super_slurm_%u_optim_h.sh',i),'w');
        vcount = 0;
        for vb=1:VBMAX

            outfile_fullname = sprintf('%s_vbatch%uof%u_Output_Roi_byparc.mat',outfile,vb,VBMAX);
            outfile_consname = sprintf('%s_Output_Roi_byparc.mat.mat',outfile);

            if ~exist(outfile_fullname,'file') && ~exist(outfile_consname,'file')

            vcount = vcount+1;

            filename = sprintf('tex_sbatch_s%u_%u_optim_h.sl',i,vb);
            fid = fopen(filename,'w');
            % header prep stuff
            fprintf(fid,'#!/bin/bash -l\n');
            fprintf(fid,'#SBATCH --job-name=tex_sbatch_s%u_%u_optim\n',i,vb);
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
            % job-specific roi subset
            roiixes = fullixes( (5*(vb-1)+1):min([5*vb numel(fullixes)]) );
            roiixlist=sprintf('%u,',roiixes);
            roiixlist=roiixlist(1:end-1);
            % write job command
            fprintf(fid,'matlab -nodisplay -nojvm -singleCompThread -r "addpath /lustre03/project/6071833/nwc/NIFTI_tools; addpath /lustre03/project/6071833/nwc/tex_code; texmap_wrapper3( ''%s'', ''%s'', [%s], 0, ''%s_vbatch%uof%u'', [], [], [], [], ''HDE'', ''optim'' ); exit;"\n',...
                t1_file,mskfile,roiixlist,outfile,vb,VBMAX);
            % close out the file
            fclose(fid);

            fprintf(fido,'sbatch %s\n',filename);
            else
                fprintf('skipping "%s", file already exists\n\n',outfile_fullname);
            end
        end
        fclose(fido);
        if vcount==0
            disp('no missing files - removing jobfile')
            unix(sprintf('rm super_slurm_%u_optim.sh',i))
        else
            fprintf('found %u missing files! to run!',vcount)
        end
    else
        disp('not on list!');
    end
end
disp('donne');
