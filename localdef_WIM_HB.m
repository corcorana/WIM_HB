%% set local paths to raw data & MATLAB toolboxes

% Melbourne dataset [download from https://osf.io/ey3ca/]
path_mbi_raweeg = 'path_to_Raw_EEG_Data';
path_mbi_raweye = 'path_to_Raw_Eye-Tracking_Data';
path_mbi_behav  = 'path_to_Behavioural_Data';

% Paris dataset [download from ???]
path_pbi_rawdat = 'path_to_Paris_Data';

% MATLAB path & toolboxes
path_mat    = 'C:\Users\corcoraa\Documents\MATLAB\';
path_eeglab = [path_mat, 'eeglab2024.2.1']; 		% https://sccn.ucsd.edu/eeglab/
path_heplab = [path_mat, 'HEPLAB-1.0.1']; 		% https://pandelisperakakis.info/heplab/
path_lmeEEG = [path_mat, 'lmeEEG-main']; 		% https://github.com/antovis86/lmeEEG
path_gcmi   = [path_mat, 'gcmi-0.4']; 			% https://github.com/robince/gcmi
path_edf2m  = [path_mat, 'edf-converter-1.21.0']; 	% https://github.com/uzh/edf-converter
path_TFCE   = [path_mat, 'ept_TFCE-matlab-master']; 	% https://github.com/Mensen/ept_TFCE-matlab
path_nmd    = [path_mat, 'NMDv2-00']; 			% http://www.physics.lancs.ac.uk/research/nbmphysics/diats/nmd/
path_pncst  = [path_mat, 'PhysioNet-Cardiovascular-Signal-Toolbox']; 	% https://physionet.org/content/pcst/1.0.0/

% Directory for preprocessed files
wim_preproc = fullfile(pwd, 'preproc_dat');
if ~exist(wim_preproc, 'dir')
    warning(['Unable to find ', wim_preproc, '... creating new directory'])
    mkdir(wim_preproc)
end

