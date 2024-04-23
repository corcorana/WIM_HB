%%
clear all
close all

% minimal preprocessing of raw EEG/ECG data
run WIM_HB_preproc_min

% EEG preprocessing
run WIM_HB_preproc_CLN
run WIM_HB_preproc_ASR
run WIM_HB_preproc_ICA
run WIM_HB_preproc_ICA_rej

% ECG preprocessing
run WIM_HB_preproc_IBI

