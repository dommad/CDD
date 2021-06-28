cdd_params_Elite.yaml - parameters of CDD for each charge state used in CDD-based FDR control

comet_t - Comet parameters applied in target database searches

comet_td - Comet parameters applied in target-decoy database searches

decoy_original.py - processing results of CDD- and PeptideProphet-based analysis, outputting final validation results 

execute_all.bash - the main Bash script executing all necessary codes of validation study

fragger.params - MSFragger parameters file used during construction of ground truth

pools_pep.csv - list of synthetic peptides in ProteomeTools (PT) data set used

refine_comet_msfrag_natural.py - script comparing results of Comet and MSFragger searches for HeLa data set, generating IDs of spectra to be modified (negative samples) and unchanged (positive samples)

refine_comet_msfrag_synthetic.py - script comparing results of Comet and MSFragger searches for PT data set, generating IDs of spectra to be modified (negative samples) and unchanged (positive samples)

val_files_pmz.py - script generating ground truth files to be used in the validation procedure
