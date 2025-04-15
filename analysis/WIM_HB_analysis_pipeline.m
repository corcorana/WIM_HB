%%
clear
close

% IBI, pupil, & behav measures (for stats & figs, run corresponding R script)
run WIM_HB_IBI_pup_behav

% HEP analysis
run WIM_HB_HEP
run WIM_HB_HEP_lmeEEG_state
run WIM_HB_HEP_lmeEEG_surrog

% GCMI analysis
run WIM_HB_GCMI_NMD
run WIM_HB_GCMI_lmeEEG

