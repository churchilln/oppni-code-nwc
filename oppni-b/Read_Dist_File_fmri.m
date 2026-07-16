function DistStruct = Read_Dist_File_fmri( distfile )
%
% script to read distortion-correction parameter file
% interprets bracketed KEY=[VALUE] fields for compatibility with OPPNI code
%

DistStruct=[];

if( ~isempty(distfile) )

    fid = fopen(distfile);
    if fid==-1
        error('cannot open dist file:\n\t%s\n',distfile);
    end
    tline = fgetl(fid);
    if ~ischar(tline)
        error('dist file is empty:\n\t%s\n',distfile);
    end

    %% read input file
    while ischar(tline)

        % [] bracketing fields
        istart = strfind(tline,'[')+1;
        iend   = strfind(tline,']')-1;

        if( numel(istart)==1 && numel(iend)==1 )

            ieq = strfind(tline,'=');
            if isempty(ieq)
                fieldname = upper(strtrim(tline(1:(istart-2))));
            else
                fieldname = upper(strtrim(tline(1:(ieq(1)-1))));
            end
            fieldname(fieldname==' ') = [];
            fieldval = strtrim(tline(istart:iend));

            if ~isempty(fieldname)
                DistStruct.(fieldname) = fieldval;
            end
        end

        tline = fgetl(fid);
        if isempty(tline)
            tline = fgetl(fid);
        end
    end
    fclose(fid);

end
