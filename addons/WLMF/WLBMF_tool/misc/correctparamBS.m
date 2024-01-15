function [paramBS] = correctparamBS(varargin)
% function [paramBS] = correctparamBS(varargin)
% write correct / non redundant bootstrap parameters for BSresample and
% ndimBSresample
% 2 possible usages:
%% A) [paramBS] = correctparamBS(paramBS)
%       checks bootstrap parameters and corrects them, if necessary
%% B) [paramBS] = correctparamBS(B1, B2, Blocklength, Method)
%       writes correct bootstrap parameters
%
% Herwig Wendt, Lyon, 2006 - 2008

%% STEP 1: READ PARAMETERS
if nargin==1
    STRUCT=1;
    tempparamBS=varargin{1}; 
    InputError1=('Error: single input argument must be structure with 9 fields: B1, B2, Blocklength, Method');
    InputError2=('Error: paramBS must contain 4 fields: B1, B2, Blocklength, Method');
    try  fnames=fieldnames(tempparamBS); catch error(InputError1, 0); end
    try  B1=getfield(tempparamBS, fnames{1}); catch error(InputError2, 0); end
    try  B2=getfield(tempparamBS, fnames{2}); catch error(InputError2, 0); end
    try  BlockLength=getfield(tempparamBS, fnames{3}); catch error(InputError2, 0); end
    try  METHOD=getfield(tempparamBS, fnames{4}); catch error(InputError2, 0); end    
elseif nargin==4
    STRUCT=0;
    B1=varargin{1};
    B2=varargin{2};
    BlockLength=varargin{3};
    METHOD=varargin{4};
else
    error('Wrong number of input arguments');
end

Method=[];
if find(METHOD==1); Method=[Method 1]; end
if find(METHOD==2); Method=[Method 2]; end
if find(METHOD==3); Method=[Method 3]; end
if find(METHOD==4); Method=[Method 4]; end
if find(METHOD==5); Method=[Method 5]; end
if find(METHOD==6); Method=[Method 6]; end

%% STEP 2: CHECK PARAMETERS
% block length
if BlockLength<1; BlockLength=1; end

paramBS=struct('Nresamp1',B1,'Nresamp2',B2, 'blocklength', BlockLength, 'Method', Method);

% check if double bootstrap methods can be calculated
if ~(B2>1) 
    Method=Method(find(Method~=4));
    Method=Method(find(Method~=5));
    Method=Method(find(Method~=6));
    paramBS=struct('Nresamp1',B1,'Nresamp2',B2, 'blocklength', BlockLength, 'Method', Method);
end

% Check if double bootstrap necessary
if ~(length(find(Method==4))|length(find(Method==5))|length(find(Method==6))); 
    B2=1; 
    paramBS=struct('Nresamp1',B1,'Nresamp2',B2, 'blocklength', BlockLength, 'Method', Method);
end