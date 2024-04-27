%% set local paths to raw data & MATLAB toolboxes

% Behavioural & raw EEG data -- download from https://osf.io/ey3ca/
path_raweeg = 'path_to_Raw_EEG_data';
path_behav  = 'path_to_Behavioural_data';

% MATLAB path & toolboxes (update versions as required)
path_mat    = 'path_to_MATLAB';
path_eeglab = [path_mat, 'eeglab2019.1'];
path_dipfit = dir(fullfile(path_eeglab, 'plugins', 'dipfit*'));
path_pncst  = [path_mat, 'PhysioNet-Cardiovascular-Signal-Toolbox'];
path_heplab = [path_mat, 'HEPLAB-1.0.1'];
path_gcmi   = [path_mat, 'gcmi-0.4'];
path_lmeEEG = [path_mat, 'lmeEEG-main'];
path_TFCE   = [path_mat, 'ept_TFCE-matlab-master'];
path_ftrip  = [path_mat, 'fieldtrip-20220810'];

% derivatives
eeg_preproc = ['preproc' filesep 'preprocEEG'];
if ~exist(eeg_preproc, 'dir')
    warning(['Unable to find ', eeg_preproc, '... creating new directory'])
    mkdir(eeg_preproc)
end

