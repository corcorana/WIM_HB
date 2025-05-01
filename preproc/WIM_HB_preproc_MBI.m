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

raw_eeg = dir([path_mbi_raweeg filesep '*.vhdr']);

%% batch process

subjs = readmatrix(['..' filesep 'subjs.txt']);

for ix = 1: length(subjs)

    snum = subjs(ix);
    sname =  ['MWI', num2str(snum)];
    fprintf('Loading %s...\n', sname)

    % load raw EEG / ECG
    fname = raw_eeg( contains({raw_eeg.name}, sname) ).name;
    EEG = pop_loadbv( path_mbi_raweeg, fname, [], [1:63, 65] );
    
    % load behav log
    behav_name = dir([path_mbi_behav filesep '*s' num2str(snum) '*.mat']);
    load([path_mbi_behav filesep behav_name.name]);

    % append vars to EEG struct
    EEG.probe_res = probe_res;
    EEG.test_res = test_res;

    % import channel locations and specify online ref (AFz)
    EEG = pop_chanedit(EEG, 'append', EEG.nbchan, 'changefield', {EEG.nbchan+1, 'labels', 'AFz'},...
        'lookup', fullfile( path_eeglab, 'plugins', 'dipfit', 'standard_BESA', 'standard-10-5-cap385.elp'),...
        'changefield', {EEG.nbchan, 'type', 'ECG'}, 'setref', {1:EEG.nbchan-1, 'AFz'} );

    % compute average ref & add online ref back into array
    EEG = pop_reref( EEG, [], 'refloc', struct('labels', {'AFz'},'type',{'EEG'},...
       'theta',{0},'radius',{0.37994},'X',{79.0255},'Y',{0},'Z',{31.3044},...
       'sph_theta',{0},'sph_phi',{21.61},'sph_radius',{85},...
       'urchan',{EEG.nbchan+1},'ref',{''},'datachan',{0}), 'exclude', find(strcmp({EEG.chanlocs.type}, 'ECG')) );
  
    % reref to TP9/10
    EEG = pop_reref( EEG, [find(strcmp({EEG.chanlocs.labels}, 'TP9')), find(strcmp({EEG.chanlocs.labels}, 'TP10'))],...
        'exclude', find(strcmp({EEG.chanlocs.type}, 'ECG')) ); 
    
    % trim off edges
    bsta = find(ismember({EEG.event.type}, 'B  1'));
    bsta = bsta(end-5);
    bend = find(ismember({EEG.event.type}, 'K  1'));
    bend = bend(end);
    EEG = pop_select( EEG, 'point', [EEG.event(bsta).latency-EEG.srate*6, EEG.event(bend).latency+EEG.srate*10] );
    
    % save .set files
    EEG = pop_saveset( EEG, 'filename', [sname, '_min.set'], 'filepath', proc_min);

end
