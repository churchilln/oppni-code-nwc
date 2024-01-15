function File_Existence_Checker_perf( InputStruct,outpath,stage )
% *PERF
% File_Existence_Checker_perf( InputStruct,stage ) --> verifies if appropriate data exists, allowing us to proceed
%                                                 compatible with gzipped and unzipped formats 
%

if stage==0

    for ns=1:numel(InputStruct)

        lstr = sprintf('Missing File (line %u/%u)! ', ns,numel(InputStruct));

        % perf n physio n task
        for nr=1:InputStruct(ns).N_perf
            lstr = sprintf('Missing File (line %u/%u, perf-run %u/%u)! ', ns,numel(InputStruct),nr,InputStruct(ns).N_perf);
            if ~exist( InputStruct(ns).prun(nr).PERF_filename ,'file')
                error('%sExpected to find PERF file but did not:\n\t%s\n', lstr, InputStruct(ns).prun(nr).PERF_filename );
            end
            if ~isempty(InputStruct(ns).prun(nr).TASK_filename) && ~exist( InputStruct(ns).prun(nr).TASK_filename ,'file')
                error('%sExpected to find TASK file but did not:\n\t%s\n', lstr, InputStruct(ns).prun(nr).TASK_filename );
            end
        end
        for nr=1:InputStruct(ns).N_m0ref
            lstr = sprintf('Missing File (line %u/%u, m0ref-run %u/%u)! ', ns,numel(InputStruct),nr,InputStruct(ns).N_perf);
            if strcmpi(InputStruct(ns).m0loc,'separate') && ~exist(InputStruct(ns).mrun(nr).M0REF_filename ,'file')
                error('%sExpected to find M0REF file but did not:\n\t%s\n', lstr, InputStruct(ns).prun(nr).M0REF_filename );
            end
        end
        % struct
        nr=1;
        lstr = sprintf('Missing File (line %u/%u, anat-run %u/%u)! ', ns,numel(InputStruct),nr,InputStruct(ns).N_anat);
        if ~exist( InputStruct(ns).arun(nr).ANAT_filename ,'file')
            error('%sExpected to find ANAT file but did not:\n\t%s\n', lstr, InputStruct(ns).arun(nr).ANAT_filename );
        end
    end

elseif stage==1

    for ns=1:numel(InputStruct)
        
        opath0 = [outpath,'/',InputStruct(ns).PREFIX,'/rawdata'];
        if ~exist(opath0,'dir')
            error('Missing Directory (line %u/%u)! Expected to find raw data folder but did not:\n\t%s\n',opath0)
        end
        
        for nr=1:InputStruct(ns).N_perf

            lstr = sprintf('Missing File (line %u/%u, perf-run %u/%u)! ', ns,numel(InputStruct),nr,InputStruct(ns).N_perf);
            if ~exist( sprintf('%s/perf%u.nii',opath0,nr), 'file') && ~exist( sprintf('%s/perf%u.nii.gz',opath0,nr), 'file')
                error('%sExpected to find raw data perf file but did not:\n\t%s\n', lstr, sprintf('%s/perf%u.nii(.gz)',opath0,nr) );
            end
            if ~strcmpi( InputStruct(ns).PWDROP, 'NONE' ) && ~exist( sprintf('%s/perf%u_drop.nii',opath0,nr), 'file') && ~exist( sprintf('%s/perf%u_drop.nii.gz',opath0,nr), 'file')
                error('%sExpected to find raw data perf file but did not:\n\t%s\n', lstr, sprintf('%s/perf%u_drop.nii(.gz)',opath0,nr) );
            end
        end
        if ~exist( sprintf('%s/m0ref_cat.nii',opath0), 'file') && ~exist( sprintf('%s/m0ref_cat.nii.gz',opath0), 'file')
            error('%sExpected to find raw data m0ref file but did not:\n\t%s\n', lstr,sprintf('%s/m0ref_cat.nii(.gz)',opath0) );
        end

        % port over anatomical data
        nr=1;
        lstr = sprintf('Missing File (line %u/%u, anat-run %u/%u)! ', ns,numel(InputStruct),nr,InputStruct(ns).N_anat);
        if ~exist(sprintf('%s/anat%u.nii',opath0,nr),'file') && ~exist(sprintf('%s/anat%u.nii.gz',opath0,nr),'file')
            error('%sExpected to find raw data anat file but did not:\n\t%s\n', lstr, sprintf('%s/anat%u.nii(.gz)',opath0,nr) );
        end
        if ( (ischar(InputStruct(ns).arun(nr).ZCLIP_thr) && strcmpi(InputStruct(ns).arun(nr).ZCLIP_thr,'AUTO')) || ... % either auto-clip or specified numeric value
           (isnumeric(InputStruct(ns).arun(nr).ZCLIP_thr) && isfinite(InputStruct(ns).arun(nr).ZCLIP_thr)) ) && ...
            ~exist(sprintf('%s/anat%u_zclip.nii',opath0,nr),'file') && ...
            ~exist(sprintf('%s/anat%u_zclip.nii.gz',opath0,nr),'file')
            error('%sExpected to find raw data anat file but did not:\n\t%s\n', lstr, sprintf('%s/anat%u_zclip.nii(.gz)',opath0,nr) );
        end
    end
end

