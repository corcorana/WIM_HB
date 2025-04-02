%%
clear
close

% minimal preprocessing of raw EEG/ECG data from both sites
run WIM_HB_preproc_MBI
run WIM_HB_preproc_PBI

% EEG-ECG artefact rejection
run WIM_HB_preproc_CLN
run WIM_HB_preproc_ASR
run WIM_HB_preproc_ECG
run WIM_HB_preproc_ICA
run WIM_HB_preproc_ICA_rej
run WIM_HB_preproc_CFA
run WIM_HB_preproc_BHI

% EYE preprocessing
run WIM_HB_preproc_EYE
