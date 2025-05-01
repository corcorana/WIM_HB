%%
clear
close

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))

proc_min = fullfile(wim_preproc, 'min');
if ~exist(proc_min, 'dir')
    mkdir(proc_min)
end

raw_eeg = dir(['D:\WIM\PBI\sub*' filesep '*.vhdr']);
behav = dir(['D:\WIM\PBI\sub*' filesep '*behavres*.mat']);

%% batch process

for ix = 1: length(raw_eeg)

    snum = ix;
    sname = sprintf('MWI%03d', snum);
    fprintf('Loading %s...\n', sname)

    % load raw EEG
    if ix == 2
        fname = raw_eeg( contains({raw_eeg.name}, '07-12_MP' ) ).name;
    else
        fname = raw_eeg( contains({raw_eeg.name}, sprintf('_%03d', snum) ) ).name;
    end
    EEG = pop_loadbv( raw_eeg(ix).folder, fname, [], [1:63, 67] );
 
    % downsample to match MBI data
    EEG = pop_resample(EEG, 500);

    % load behav log
    if ix == 18    
        behav_name = behav( contains({behav.name}, 'sMR' )).name;
    else
        behav_name = behav( contains({behav.name}, sprintf('s%03d', snum) )).name;
    end
    load([behav(ix).folder filesep behav_name]);

    % append vars to EEG struct
    EEG.probe_res = probe_res;
    EEG.test_res = test_res;

    % import channel locations and specify online ref (FCz)
    EEG = pop_chanedit(EEG, 'append', EEG.nbchan, 'changefield', {EEG.nbchan+1, 'labels', 'FCz'},...
        'lookup', fullfile( path_eeglab, 'plugins', 'dipfit', 'standard_BESA', 'standard-10-5-cap385.elp'),...
        'changefield', {EEG.nbchan, 'type', 'ECG'}, 'setref', {1:EEG.nbchan-1, 'FCz'} );

    % compute average ref & add online ref back into array
    EEG = pop_reref( EEG, [], 'refloc', struct('labels', {'FCz'},'type',{'EEG'},...
       'theta',{0},'radius',{0.12662},'X',{32.9279},'Y',{0},'Z',{78.363},...
       'sph_theta',{0},'sph_phi',{67.208},'sph_radius',{85},'urchan',{EEG.nbchan+1},'ref',{''},'datachan',{0}),...
       'exclude', find(strcmp({EEG.chanlocs.type}, 'ECG')) );
  
    % reref to TP9/10
    EEG = pop_reref( EEG, [find(strcmp({EEG.chanlocs.labels}, 'TP9')), find(strcmp({EEG.chanlocs.labels}, 'TP10'))],...
        'exclude', find(strcmp({EEG.chanlocs.type}, 'ECG')) ); 
 
    % trim off edges
    f = find(strcmp({EEG.event.type}, 'S  2'));
    d = diff(find(strcmp({EEG.event.type}, 'S  2')))==1;
    bsta = f(d);
    f = find(strcmp({EEG.event.type}, 'S 22'));
    bend = f(f>bsta(1));
    EEG = pop_select( EEG, 'point', [EEG.event(bsta(1)).latency-EEG.srate*6, EEG.event(bend(end)).latency+EEG.srate*10] );

    % strip original filename
    EEG.comments = [];

    % save .set files
    EEG = pop_saveset( EEG, 'filename', [sname, '_min.set'], 'filepath', proc_min);

end
