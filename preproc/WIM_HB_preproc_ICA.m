%%
clear all
close all

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))

proc_ica = fullfile(eeg_preproc, 'ica');
if ~exist(proc_ica, 'dir')
    mkdir(proc_ica)
end

%% batch process

subjs = dlmread(['..' filesep 'subjs.txt']);

for ix = 1:length(subjs)
    snum = subjs(ix);
    sname = ['MWI' num2str(snum)];
    fprintf(['Loading ' sname '...\n'])

    % load ASR'd data
    EEG = pop_loadset('filename', [sname, '_asr.set'], 'filepath', fullfile(eeg_preproc, 'asr'));

    % drop ECG
    EEG = pop_select( EEG, 'nochannel', {'ECG'});

    % low-pass filter
    EEG = pop_eegfiltnew(EEG, [], 40);
    
    % interpolate missing channels
    EEG = pop_interp(EEG, EEG.urchanlocs, 'spherical');
    
    % average reference
    EEG = pop_reref( EEG, [],'refloc',struct('labels',{'TP9','TP10'},'type',{'EEG','EEG'},...
        'theta',{-108.393,108.393},'radius',{0.66489,0.66489},'X',{-23.3016,-23.3016},'Y',{70.0758,-70.0758},...
        'Z',{-42.0882,-42.0882},'sph_theta',{108.393,-108.393},'sph_phi',{-29.68,-29.68},'sph_radius',{85,85},...
        'urchan',{10,21},'ref',{'TP9 TP10','TP9 TP10'},'datachan',{0,0}));

    % run ICA
    EEG = pop_runica(EEG, 'icatype', 'binica', 'extended',1, 'interrupt','off');
    
    % class ICs
    EEG = iclabel(EEG, 'default');
    pop_viewprops(EEG, 0, 1:35, {'freqrange', [1 45]}); % view component properties

    % save decomposition
    EEG = pop_saveset( EEG, 'filename', [sname, '_ica.set'], 'filepath', proc_ica );
    saveas(gcf, fullfile(eeg_preproc, 'ica', [sname, '_ica_topo.png']) );
    close
    
end

% clear temporary files
delete binica*