%%
clear
close

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))

%% batch process

subjs = dir(fullfile(wim_preproc, 'ica', 'MWI*_ica.set'));

autoICA = 1;    % perform IC rejection automatically with ICLabel (1 = yes, else manual review)
crit = .9;      % threshold for ICLabel rejection
rejClass = {'Eye', 'Heart', 'Channel Noise'}; % component classes to reject

fname = sprintf('WIM_HB_ICA_rejComps_%s.csv', datestr(now, 'yyyy-mm-dd'));
fid = fopen( fname, 'w' );
fprintf( fid, '%s, %s, %s, %s\n', 'Subject', rejClass{:});
fclose( fid );
rej_comps = cell(length(subjs), 3);

for ix = 1:length(subjs)

    sname = subjs(ix).name(1:6);
    fprintf(['Loading ' sname '...\n'])

    % load ICA'd data
    EEG = pop_loadset('filename', [sname, '_ica.set'], 'filepath', fullfile(wim_preproc, 'ica'));
    rej_comps{ix,1} = sname; 

    % evaluate components using ICLabel
    if autoICA == 1
        lab = matches(EEG.etc.ic_classification.ICLabel.classes, rejClass);
        EEG.reject.gcompreject = sum(EEG.etc.ic_classification.ICLabel.classifications(:,lab) > crit, 2);
        rej_comps{ix,3} = sum(EEG.etc.ic_classification.ICLabel.classifications(:,lab) > crit);
        comprej = find(EEG.reject.gcompreject);
    else
        pop_viewprops(EEG, 0, 1:35, {'freqrange', [1 45]}); % view component properties
        pause
        comprej = inputdlg('Components for rejection:','',[1 35]);
        close
        comprej = str2num(comprej{1});
    end
    
    rej_comps{ix,2} = comprej;
    EEG = pop_subcomp(EEG, comprej);

    % save ICA-corrected data
    EEG = pop_saveset( EEG, 'filename', [sname, '_ica_rej.set'], 'filepath', fullfile(wim_preproc, 'ica'));

    % record number of rejected components per subject
    fid = fopen( fname, 'a' );
    fprintf( fid, '%s, %d, %d, %d\n', sname, rej_comps{ix,3} );
    fclose( fid );
    save(sprintf('WIM_HB_ICA_rejComps_%s.mat', datestr(now, 'yyyy-mm-dd')), 'rej_comps')

end
