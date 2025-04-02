%%
clear
close

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))

proc_ica = fullfile(wim_preproc, 'ica');
if ~exist(proc_ica, 'dir')
    mkdir(proc_ica)
end

%% batch process

subjs = dir(fullfile(wim_preproc, 'ecg', '*_qrs_screened.set'));

for ix = 1 :length(subjs)

    sname = subjs(ix).name(1:6);
    fprintf(['Loading ' sname '...\n'])

    % load preprocessed EEG & ECG data
    EEG = pop_loadset('filename', [sname, '_asr.set'], 'filepath', fullfile(wim_preproc, 'asr'));
    ECG = pop_loadset('filename', subjs(ix).name, 'filepath', fullfile(wim_preproc, 'ecg') );

    % combine data
    EEG.event = ECG.event;
    EEG.data = [EEG.data; ECG.data];
    EEG.nbchan = EEG.nbchan+1;
    EEG.chanlocs = [EEG.chanlocs, ECG.chanlocs];
    
    % downsample 
    EEG = pop_resample(EEG, 250);

    % run ICA
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1, 'interrupt','off');

    % class ICs
    EEG = iclabel(EEG, 'default');
    
    % visualise & save topos
    pop_viewprops(EEG, 0, 1:35, {'freqrange', [1 45]}); % view component properties
    saveas(gcf, fullfile(wim_preproc, 'ica', [sname, '_ica_topo.png']) );
    close
   
    % save decomposition
    EEG = pop_saveset( EEG, 'filename', [sname, '_ica.set'], 'filepath', proc_ica );

end

