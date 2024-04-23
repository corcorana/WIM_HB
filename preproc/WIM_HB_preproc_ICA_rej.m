%%
clear all
close all

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))

%% batch process

subjs = dlmread(['..' filesep 'subjs.txt']);

bad_comps = cell(length(subjs), 3);

autoICA = 1;    % perform IC rejection automatically with ICLabel (1 = yes, else manual review)
crit = .9;      % threshold for ICLabel rejection

for ix = 1:length(subjs)
    snum = subjs(ix);
    sname = ['MWI' num2str(snum)];
    fprintf(['Loading ' sname '...\n'])

    % load ICA'd data
    EEG = pop_loadset('filename', [sname, '_ica.set'], 'filepath', fullfile(eeg_preproc, 'ica'));

    % estimate amount of variance explained by each IC timeseries
    compseries = eeg_getdatact(EEG, 'component', 1:size(EEG.icaweights,1));
    pvar = var(compseries,[],2)/sum(var(compseries,[],2))*100;
    bad_comps{ix,1} = sum(pvar>1);
         
    % evaluate components using ICLabel
    if autoICA == 1
        lab = matches(EEG.etc.ic_classification.ICLabel.classes, {'Eye', 'Heart', 'Channel Noise'});
        EEG.reject.gcompreject = sum(EEG.etc.ic_classification.ICLabel.classifications(pvar>1,lab) > crit, 2);
        bad_comps{ix,3} = sum(EEG.etc.ic_classification.ICLabel.classifications(pvar>1,lab) > crit);
        comprej = find(EEG.reject.gcompreject);
    else
        pop_viewprops(EEG, 0, 1:28, {'freqrange', [1 45]}); % view component properties
        pause
        comprej = inputdlg('Components for rejection:','',[1 35]);
        close
        comprej = str2num(comprej{1});
    end
    bad_comps{ix,2} = comprej;
    EEG = pop_subcomp(EEG, comprej);
    fprintf(['...Rejecting ' num2str(length(comprej)) ' components...\n'])

    % save ICA-corrected data
    EEG = pop_saveset( EEG, 'filename', [sname, '_ica_rej.set'], 'filepath', fullfile(eeg_preproc, 'ica'));
    save('WIM-HB_ICA_bad_comps.mat', 'bad_comps')
end

% save summary of bad_comps
ncomps = cellfun('length', bad_comps);
dlmwrite('WIM-HB_ICA_bad_comps.txt', [subjs, ncomps(:,2)] )
