%%
clear
close

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))

proc_bhi = fullfile(wim_preproc, 'bhi');
if ~exist(proc_bhi, 'dir')
    mkdir(proc_bhi)
end

%% batch process
subjs = dir(fullfile(wim_preproc, 'ica', 'MWI*_ica_rej_cfa.set'));

for ix = 1 :length(subjs)
    
    sname = subjs(ix).name(1:6);
    fprintf(['Loading ' sname '...\n'])

    % load CFA-corrected EEG 
    EEG = pop_loadset( 'filename', [sname, '_ica_rej_cfa.set'], 'filepath', fullfile(wim_preproc, 'ica'));
    fs = EEG.srate;

    % remove CFA-corrected ECG
    EEG = pop_select( EEG, 'rmchantype', {'ECG'});

    % load intact ECG
    ECG = pop_loadset( 'filename', [sname, '_ica.set'], 'filepath', fullfile(wim_preproc, 'ica'));
    ECG = pop_select( ECG, 'chantype', {'ECG'});

    % insert NaNs in place of missing channel data
    miss_chan = setdiff({EEG.urchanlocs.labels}, {EEG.chanlocs.labels});
    if ~isempty(miss_chan)
        EEG = pop_interp(EEG, EEG.urchanlocs, 'spherical');
        EEG.data(matches({EEG.chanlocs.labels}, miss_chan), :) = nan( numel(miss_chan), EEG.pnts );
    end

    % ensure consistent channel order across datasets
    [lab,idx] = sort({EEG.chanlocs.labels});
    EEG.chanlocs = EEG.chanlocs(idx);
    EEG.data = EEG.data(idx,:);

    % reunite with ECG
    EEG.data = [EEG.data; ECG.data];
    EEG.chanlocs = [EEG.chanlocs, ECG.chanlocs];
    EEG.nbchan = numel(EEG.chanlocs);

    % save 
    EEG = pop_saveset( EEG, 'filename', [sname, '_bhi.set'], 'filepath', proc_bhi );

end

