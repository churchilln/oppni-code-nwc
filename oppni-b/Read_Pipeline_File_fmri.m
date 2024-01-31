function PipeStruct = Read_Pipeline_File_fmri( pipelinefile )
%
% PipeStruct.(field) = cell array of strings (allows multiple pipeline optionts to test)

PipeStruct=[];
if ~contains(pipelinefile,'=')
    fid   = fopen(pipelinefile);
    if fid==-1
        error('cannot open pipeline file:\n\t%s\n',pipelinefile);
    end
    newline = fgetl(fid);
    if ~ischar(newline)
        error('pipeline file is empty:\n\t%s\n',pipelinefile);
    end
    newline(newline==' ') = [];
    pipelinestring=[];
    while ischar(newline)
        %
        pipelinestring = [pipelinestring ' ' newline];
        newline      = fgetl(fid);
        newline(newline==' ') = [];
    end
    fclose(fid);
else
    pipelinestring = pipelinefile; %% just directly assign the string
end

steplist = {'PNAME','AMASK','AWARP','ASEG','DESPIKE','RICOR','TSHIFT','FWARP','SMOOTH','DETREG','GSREG','MOTREG','ROIREG','TASKREG','LOPASS','COMPFILT'}; % ,'PHYCOR'

ileft  = strfind( pipelinestring, '[' );
iright = strfind( pipelinestring, ']' );
for(s=1:length(steplist))
    iStep = strfind( upper(pipelinestring), strcat(steplist{s},'=') );
    if isempty(iStep)
        if strcmpi(steplist,'PNAME')
            error('Need to give this pipeline an identifier "PNAME". This is a unique identifier for your pipeline of interest.')
        else
            % for saftety, we currently force the user to explicitly define what they want to do with the pipeline step
            error('Pipeline step %s has not been defined, please check the spelling!\n',steplist{s});
            %PipeStruct.(steplist{s}) = {'0'};
        end
    else
        bleft       = ileft (ileft >iStep);   bleft= bleft(1);
        bright      = iright(iright>iStep);  bright=bright(1);
        pipeargs    = regexp(  pipelinestring((bleft+1):(bright-1))  ,',','split');
        if s>1 && numel(pipeargs)==1 % second field added-
            pipeargs{2} = [];
        end
        if strcmpi(steplist{s},'PNAME') && numel(pipeargs)>1
            error('Pipeline step %s has multiple comma-separated arguments\n only allowed to specify one argument!',steplist{s});
        end
        if strcmpi(pipeargs{1},'0') % alternate coding for OFF
            pipeargs{1} = 'OFF';
        end
        PipeStruct.(steplist{s}) = pipeargs;
    end
end

% subpipe - anatwarp
a=sprintf('%s_',PipeStruct.AMASK{:});
b=sprintf('%s_',PipeStruct.AWARP{:});
    tmp_Warp_ID = ['#Warp#','.',a(1:end-1),'.',b(1:end-1)];
% subpipe - segmentation
a=sprintf('%s_',PipeStruct.ASEG{:});
    tmp_Seg_ID = ['#Seg#','.',a(1:end-1)];
% subpipe - func-pre+reg
a=sprintf('%s_',PipeStruct.DESPIKE{:});
b=sprintf('%s_',PipeStruct.RICOR{:});
c=sprintf('%s_',PipeStruct.TSHIFT{:});
d=sprintf('%s_',PipeStruct.FWARP{:});
e=sprintf('%s_',PipeStruct.SMOOTH{:});
    tmp_P1_ID = ['#P1#','.',a(1:end-1),'.',b(1:end-1),'.',c(1:end-1),'.',d(1:end-1),'.',e(1:end-1)];
% subpipe - func-everything else
a=sprintf('%s_',PipeStruct.DETREG{:});
b=sprintf('%s_',PipeStruct.GSREG{:});
c=sprintf('%s_',PipeStruct.MOTREG{:});
d=sprintf('%s_',PipeStruct.ROIREG{:});
e=sprintf('%s_',PipeStruct.TASKREG{:});
f=sprintf('%s_',PipeStruct.LOPASS{:});
g=sprintf('%s_',PipeStruct.COMPFILT{:});
    tmp_P2_ID = ['#P2#','.',a(1:end-1),'.',b(1:end-1),'.',c(1:end-1),'.',d(1:end-1),'.',e(1:end-1),'.',f(1:end-1),'.',g(1:end-1)];

%--> now rename to capture dependencies
PipeStruct.Warp_ID = [tmp_Warp_ID];
PipeStruct.Seg_ID  = [tmp_Warp_ID,'-',tmp_Seg_ID]; %- depends uniquely on mask+warp
PipeStruct.P1_ID   = [tmp_Warp_ID,'-',tmp_P1_ID]; %- depends uniquely on mask+warp too
PipeStruct.P2_ID   = [tmp_Warp_ID,'-',tmp_Seg_ID,'-',tmp_P1_ID,'-',tmp_P2_ID]; %- depends on mask+warp, segmentation, P1 steps plus P2 steps
