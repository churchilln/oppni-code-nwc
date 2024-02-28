function File_Existence_Checker_fmri_noproc( InputStruct,outpath,stage )
%
% File_Existence_Checker_fmri( InputStruct,stage ) --> verifies if appropriate data exists, allowing us to proceed
%                                                 compatible with gzipped and unzipped formats 
%

if stage==0

    for ns=1:numel(InputStruct)

        lstr = sprintf('Missing File (line %u/%u)! ', ns,numel(InputStruct));
        if ~isempty(InputStruct(ns).SEED_filename) && ~exist( InputStruct(ns).SEED_filename ,'file')
            error('%sExpected to find SEED file but did not:\n\t%s\n', lstr, InputStruct(ns).SEED_filename );
        end       

        % func n physio n task
        for nr=1:InputStruct(ns).N_func
    
            lstr = sprintf('Missing File (line %u/%u, func-run %u/%u)! ', ns,numel(InputStruct),nr,InputStruct(ns).N_func);
            if ~exist( InputStruct(ns).frun(nr).FUNC_filename ,'file')
                error('%sExpected to find FUNC file but did not:\n\t%s\n', lstr, InputStruct(ns).frun(nr).FUNC_filename );
            end
            if ~isempty(InputStruct(ns).frun(nr).TASK_filename) && ~exist( InputStruct(ns).frun(nr).TASK_filename ,'file')
                error('%sExpected to find TASK file but did not:\n\t%s\n', lstr, InputStruct(ns).frun(nr).TASK_filename );
            end
            if ~isempty(InputStruct(ns).frun(nr).MPE_filename) && ~exist( InputStruct(ns).frun(nr).MPE_filename ,'file')
                error('%sExpected to find MPE file but did not:\n\t%s\n', lstr, InputStruct(ns).frun(nr).MPE_filename );
            end            
            if ~isempty(InputStruct(ns).frun(nr).PHYSIO_filename)
                taglist = {'puls','resp'};
                for i=1:numel(taglist)
                    if ~exist( sprintf('%s.%s.1D',InputStruct(ns).frun(nr).PHYSIO_filename,taglist{i}) ,'file')
                        if exist( sprintf('%s.%s',InputStruct(ns).frun(nr).PHYSIO_filename,taglist{i}) ,'file')
                            error('%sFound PHYSIO %s file, but improperly formatted. Please convert to .1D format:\n\t%s\n', lstr, taglist{i}, InputStruct(ns).frun(nr).PHYSIO_filename )
                        else
                            error('%sExpected to find PHYSIO %s file but did not:\n\t%s\n', lstr, taglist{i}, InputStruct(ns).frun(nr).PHYSIO_filename )
                        end
                    end
                end
            end
        end
    end

elseif stage==1

    for ns=1:numel(InputStruct)
        
        opath0 = [outpath,'/',InputStruct(ns).PREFIX,'/rawdata'];
        if ~exist(opath0,'dir')
            error('Missing Directory (line %u/%u)! Expected to find raw data folder but did not:\n\t%s\n',opath0)
        end
        
        for nr=1:InputStruct(ns).N_func

            lstr = sprintf('Missing File (line %u/%u, func-run %u/%u)! ', ns,numel(InputStruct),nr,InputStruct(ns).N_func);
            if ~exist( sprintf('%s/func%u.nii',opath0,nr), 'file') && ~exist( sprintf('%s/func%u.nii.gz',opath0,nr), 'file')
                error('%sExpected to find raw data func file but did not:\n\t%s\n', lstr, sprintf('%s/func%u.nii(.gz)',opath0,nr) );
            end
            if ~isempty(InputStruct(ns).frun(nr).FMASK_filename) && ~exist(sprintf('%s/func%u_mask.nii.gz',opath0,nr),'file')
                error('%sExpected to find raw mask file but did not:\n\t%s\n', lstr, sprintf('%s/func%u_mask.nii.gz',opath0,nr) );
            end
            if ~isempty(InputStruct(ns).frun(nr).MPE_filename) && ~exist(sprintf('%s/func%u_mpe',opath0,nr),'file')
                error('%sExpected to find raw mpe file but did not:\n\t%s\n', lstr, sprintf('%s/func%u_mpe',opath0,nr) );
            end
            if ~isempty(InputStruct(ns).frun(nr).PHYSIO_filename) && ~exist(sprintf('%s/physio%u.resp.1D',opath0,nr),'file')
                error('%sExpected to find raw data physio file but did not:\n\t%s\n', lstr, sprintf('%s/physio%u.resp.1D',opath0,nr) );
            end
            if ~isempty(InputStruct(ns).frun(nr).PHYSIO_filename) && ~exist(sprintf('%s/physio%u.puls.1D',opath0,nr),'file')
                error('%sExpected to find raw data physio file but did not:\n\t%s\n', lstr, sprintf('%s/physio%u.puls.1D',opath0,nr) );
            end
        end
    end
end

