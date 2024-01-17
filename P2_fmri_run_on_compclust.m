function P2_fmri_run_on_compclust( inputfile, pipelinefile, paramlist, outpath, big_skip, mask_subj_idxes, subj_subset )


fido = fopen(sprintf('super_slurm_run_oppni-b.sh'),'w');

function P2_fmri_dataProcessing( inputfile, pipelinefile, paramlist, outpath, big_skip, mask_subj_idxes, subj_subset )

fprintf(fid,'matlab -nodisplay -nojvm -singleCompThread -r "addpath ??; P2_fmri_dataProcessing( %s, %s, %s, %s, 0, %s, %u ); exit;"\n',...
    inputfile, pipelinefile, paramlist, outpath, mask_subj_idxes, subnum );



% texmap_wrapper3( ''%s'', ''%s'', [%s], 0, ''%s_vbatch%uof%u'', [], [], [], [], ''HDE'', ''optim'' ); 