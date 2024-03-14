function out = LME_ph( image_mat_exp, design_mat_exp, varnames, equation )
%
% - vox x samp
% - design includes all relevant vars

% gets data dimensions
[Nvox,Ndat] = size(image_mat_exp);
 Nvar       = size(design_mat_exp);

% cats up the design names
design_names = [{'Y'},varnames];

for v=1:Nvox

    fprintf('completed %u/%u voxels...\n',v,Nvox);

    vmat = [image_mat_exp(v,:)' design_mat_exp];
    tab  = array2table( vmat, 'VariableNames',design_names);
    %---
    out_lmm  = fitlme( tab, equation ); % mixed effect model stored

    % -pull fixed effects out of model
    [~,~,c] = fixedEffects(out_lmm);
    % store effects as matrices
    %mat_effec(v,:) = c.Estimate; % coefficients of effect
    mat_tstat(v,:) = c.tStat;    % t-statistics
    mat_pvals(v,:) = c.pValue;   % p-values
end
disp('done vox!')

inam = c.Name;
% resplitting based on design
restr = equation(strfind(equation,'~')+1:end);
restr = regexp(restr,'+','split');
restr = strtrim(restr);

out.tstat   = NaN*ones(size(image_mat_exp,1),numel(restr));
out.tstat_p = NaN*ones(size(image_mat_exp,1),numel(restr));

for i=1:numel(restr)
    ix = find( strcmpi(restr{i}, c.Name) ),
    if ~isempty(ix)
        out.tstat   = mat_tstat(:,ix);
        out.tstat_p = mat_pvals(:,ix);
    end
end

out.testname = 'lme_tstat';
