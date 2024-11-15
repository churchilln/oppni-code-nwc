function ricor_AF1( Funcfile, prefix, odir, pulsfile, respfile, acqpar, ParamCell)
%
% .ricor_AF1:
% .regression of cardiac and respiratory signals' phase effects with AFNI utilities 
% .uses 3dretroicor module

if isempty(ParamCell)
    CARD_ord  = 2;
    RESP_ord  = 2;
    INTER_ord = 0;
elseif numel(ParamCell)==3
    
    ixc=find(contains(ParamCell,'CARD-'));
    ixr=find(contains(ParamCell,'RESP-'));
    ixi=find(contains(ParamCell,'INTER-'));
    
    if isempty(ixc) || isempty(ixr) || isempty(ixi)
        error('something missing in RICOR params');
    else
        CARD_ord = str2num(ParamCell{ixc}(6:end));
        RESP_ord = str2num(ParamCell{ixr}(6:end));
        INTER_ord = str2num(ParamCell{ixi}(6:end));
    end    
else
    error('need to specify 2 fields for RICOR: e.g., CARD-2,RESP-2,INTER-0')
end

if CARD_ord ~= RESP_ord
    error('card and resp orders must match for this ricor implementation');
end
if INTER_ord>0
    error('interaction mode not available for this ricor implementation');
end
MODEL_ord = CARD_ord;

pref = [odir,'/__opptmp_p2func_ricor'];

spec_case = {'alt+z','alt+z2','alt-z','alt-z2','seq+z','seq-z'};

if ~exist(sprintf('%s/%s_ricor.nii.gz',odir,prefix),'file')

    % build directory struct recursively
    unix(sprintf('mkdir -p %s',pref));

    % embed the slice timing information!
    if sum(strcmpi(acqpar.tpatt,spec_case))>0
        disp('special case arg');
        unix(sprintf('3dTcat -prefix %s/func_tinhdr.nii.gz -tpattern %s %s',pref,acqpar.tpatt,Funcfile));
    else
        disp('custom arg');
        unix(sprintf('3dTcat -prefix %s/func_tinhdr.nii.gz -tpattern @%s %s',pref,acqpar.tpatt,Funcfile));
    end
    % now ricor it
    unix(sprintf('3dretroicor -order %u -prefix %s/%s_ricor.nii.gz -card %s -resp %s %s/func_tinhdr.nii.gz',MODEL_ord,odir,prefix,pulsfile,respfile,pref));

    unix(sprintf('rm -rf %s',pref));
else
    disp('afni-ricor already exists!')
end
