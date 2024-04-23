%%
clear all
close all

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))

proc_asr = fullfile(eeg_preproc, 'asr');
if ~exist(proc_asr, 'dir')
    mkdir(proc_asr)
end

%% batch process

subjs = dlmread(['..' filesep 'subjs.txt']);

bad_chans = cell(length(subjs), 1);
for ix = 1:length(subjs)
    snum = subjs(ix);
    sname = ['MWI' num2str(snum)];
    fprintf(['Loading ' sname '...\n'])

    % load cleanlined EEG
    EEG = pop_loadset( 'filename', [sname, '_cln.set'], 'filepath', fullfile(eeg_preproc, 'cln'));

    % exclude non-trial data (between blocks & within probes)
    bsta = find(ismember({EEG.event.type}, 'B  1'));
    bsta = bsta(2:end);
    bend = find(ismember({EEG.event.type}, 'K  1'));
    bend = bend(1:end-1);
    for bx = length(bend):-1:1
        if (EEG.event(bsta(bx)).latency) - (EEG.event(bend(bx)).latency+EEG.srate*2) > EEG.srate
            EEG = pop_select( EEG, 'nopoint', [EEG.event(bend(bx)).latency+EEG.srate*2, EEG.event(bsta(bx)).latency-1] );
        end
    end
    psta = find(contains({EEG.event.type}, 'P  '));
    pend = find(matches({EEG.event.type}, 'C  1'));
    for px = length(psta):-1:1
    	EEG = pop_select( EEG, 'nopoint',[EEG.event(psta(px)).latency+EEG.srate*2, EEG.event(pend(px)).latency-1] );
    end


    % separate for ASR
    ECG = pop_select( EEG, 'channel', {'ECG'});
    ECG.chanlocs.type = 'ECG';
    ECG.chanlocs.urchan = [];
    ECG.chanlocs.ref = [];
    EEG = pop_select( EEG, 'nochannel', {'ECG'});

    
    % compute average ref & add AFz back into array
    EEG = pop_reref( EEG, [], 'refloc', struct('labels', {'AFz'},'type',{'EEG'},...
       'theta',{0},'radius',{0.37994},'X',{79.0255},'Y',{0},'Z',{31.3044},...
       'sph_theta',{0},'sph_phi',{21.61},'sph_radius',{85},'urchan',{EEG.nbchan+1},'ref',{''},'datachan',{0}));
    
    % reref to TP9/10
    EEG = pop_reref( EEG, [find(strcmp({EEG.chanlocs.labels}, 'TP9')),...
        find(strcmp({EEG.chanlocs.labels}, 'TP10'))] );
    
    % ASR to reject bad EEG channels & repair bursts
    ASR = pop_clean_rawdata(EEG, 'FlatlineCriterion',5, 'ChannelCriterion',0.8,...
        'LineNoiseCriterion','off', 'Highpass','off', 'BurstCriterion', 20,...
        'WindowCriterion','off', 'BurstRejection','off', 'Distance','Euclidian');
    
    % drop bad channels from EEG
   	bad_chans{ix,1} = setdiff({EEG.chanlocs.labels}, {ASR.chanlocs.labels});
    ASR.urchanlocs = EEG.chanlocs;

    % reunite EEG/ECG for QRS extraction and artifact rejection
    ASR.nbchan = ASR.nbchan+1;
    ASR.data = [ASR.data; ECG.data];
    ASR.chanlocs = [ASR.chanlocs, ECG.chanlocs];

    % save combined EEG/ECG dataset
    ASR = pop_saveset( ASR, 'filename', [sname, '_asr.set'], 'filepath', proc_asr);
    save('WIM-HB_ASR_bad_chans.mat', 'bad_chans')
end

% save summary of bad_chans
dlmwrite('WIM-HB_ASR_bad_chans.txt', [subjs, cellfun('length', bad_chans)] )
