function SeedStruct = Read_Seed_File_fmri( seedfile )
%
% script to read "Split-Info" file, in either .mat or .txt format
% interprets and reformats for compatibility with OPPNI code
%

SeedStruct=[];

if( ~isempty(seedfile) )

    fid = fopen(seedfile);
    if fid==-1
        error('cannot open seed file:\n\t%s\n',seedfile);
    end
    tline = fgetl(fid);
    if ~ischar(tline)
        error('seed file is empty:\n\t%s\n',seedfile);
    end
    
    n_seed     = 0;
    n_seedloc  = 0;
    n_seedtyp  = 0;
    seedlist   = {}; %% to check for duplicates

    %% read input file
    while ischar(tline) 

        % [] bracketing fields
        istart = strfind(tline,'[')+1;
        iend   = strfind(tline,']')-1;   

        if( ~isempty(istart) && ~isempty(iend) )
              
            %% seed based analysis fields...
            
            % check for condition-specific options
            if contains(upper(tline),'SEED_LABEL')
                n_seed = n_seed+1;
                if( numel(regexp( tline(istart:iend), '[^a-zA-Z0-9_]' ))>0 )
                    error('split file info: %s seed names can only include alphanumeric values and underscores.',seedfile);
                end
                SeedStruct.seed(n_seed).label = tline(istart:iend);
                seedlist = [seedlist, {SeedStruct.seed(n_seed).label}];
            end
            % 
            if contains(upper(tline),'SEED_LOCATION') % parc, comp, multi-parc, multi-comp (binary/continuous

                n_seedloc = n_seedloc+1;
                seedlocArray{n_seedloc} = tline(istart:iend);
            end      
            % 
            if contains(upper(tline),'SEED_TYPE') % || ~isempty(strfind(upper(tline),'DURATIONS')))

                n_seedtyp = n_seedtyp+1;
                seedspcArray{n_seedtyp} = tline(istart:iend);
            end                
        end

        tline = fgetl(fid);
        if isempty(tline)
            tline = fgetl(fid);
        end
    end
    fclose(fid); 

    if( length(unique(seedlist)) < length(seedlist) )
        error('duplicate conditions found in task design of %s\n',seedfile);
    end

    % seed based quality checks
    if( (n_seed == n_seedloc) && (n_seed == n_seedtyp) && (n_seedtyp == n_seedloc) )

        for(n=1:n_seed)

            SeedStruct.seed(n).location = seedlocArray{n};                
            SeedStruct.seed(n).type = seedspcArray{n};

            if ~exist(SeedStruct.seed(n).location,'file')
                error('Seed file for %s, located at %s cannot be found! Check your paths/names!', SeedStruct.seed(n).label, SeedStruct.seed(n).location);
            end
            if ~strcmpi(SeedStruct.seed(n).type,'single') && ~strcmpi(SeedStruct.seed(n).type,'multi')
                error('Seed file for %s, has invalid type %s! Should be either "single" or "multi"!', SeedStruct.seed(n).label, SeedStruct.seed(n).type);
            end
        end
    else
        error('every seed in split-info. requires a SEED_LABEL, list of SEED_LOCATION and SEED_TYPE');
    end       
end

