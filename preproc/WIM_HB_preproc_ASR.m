%%
clear
close

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))

proc_asr = fullfile(wim_preproc, 'asr');
if ~exist(proc_asr, 'dir')
    mkdir(proc_asr)
end
proc_ecg = fullfile(wim_preproc, 'ecg');
if ~exist(proc_ecg, 'dir')
    mkdir(proc_ecg)
end

%% batch process

subjs = dir(fullfile(wim_preproc, 'cln', 'MWI*.set'));

fname = sprintf('WIM_HB_ASR_badChans_%s.csv', datestr(now, 'yyyy-mm-dd'));
fid = fopen( fname, 'w' );
fprintf( fid, '%s, %s\n', 'subID', 'rejChans');
fclose( fid );
rej_chans = cell(length(subjs),1);

for ix = 1:length(subjs)
    
    sname = subjs(ix).name(1:6);
    fprintf(['Loading ' sname '...\n'])

    % load cleanlined EEG
    EEG = pop_loadset( 'filename', subjs(ix).name, 'filepath', fullfile(wim_preproc, 'cln'));

    % low-pass filter
    EEG = pop_eegfiltnew(EEG, [], 40);

    % exclude non-trial data between blocks
    if EEG.subj_info.subID(1) == '0'        % PBI data
        f = find(strcmp({EEG.event.type}, 'S  2'));
        bsta = f(1);
        EEG = pop_select( EEG, 'nopoint', [EEG.event(1).latency, EEG.event(f(1)).latency] );
        f = find(strcmp({EEG.event.type}, 'S 22'));
        bend = f(f>bsta);
        for bx = length(bend):-1:1
            try
                EEG = pop_select( EEG, 'nopoint', [EEG.event(bend(bx)).latency+EEG.srate*4, EEG.event(bend(bx)+1).latency] );
            catch
                EEG = pop_select( EEG, 'nopoint', [EEG.event(bend(bx)).latency+EEG.srate*4, size(EEG.data, 2)] );
            end
        end
    elseif EEG.subj_info.subID(1) == '3'    % MBI data
        bsta = find(ismember({EEG.event.type}, 'B  1'));
        bsta = bsta(2:end);
        bend = find(ismember({EEG.event.type}, 'K  1'));
        bend = bend(1:end-1);
        for bx = length(bend):-1:1
            if (EEG.event(bsta(bx)).latency) - (EEG.event(bend(bx)).latency+EEG.srate*2) > EEG.srate
                EEG = pop_select( EEG, 'nopoint', [EEG.event(bend(bx)).latency+EEG.srate*2, EEG.event(bsta(bx)).latency-1] );
            end
        end
    end

    if strcmp(sname, 'MWI323')
        EEG = pop_select( EEG, 'nopoint', [1041000, 1148500] ); % ? experiment paused due to technical failure
    end

    % set ECG aside (ASR still identifies as bad channel when told to ignore)
    ECG = pop_select( EEG, 'channel',{'ECG'}); 
    ECG = pop_saveset( ECG, 'filename', [sname, '_ecg.set'], 'filepath', proc_ecg);
    EEG = pop_select( EEG, 'nochannel',{'ECG'});     

    % ASR to reject bad EEG channels & repair bursts
    ASR = pop_clean_rawdata(EEG, 'FlatlineCriterion',5, 'ChannelCriterion',0.8,...
        'LineNoiseCriterion','off', 'Highpass','off', 'BurstCriterion', 20,...
        'WindowCriterion','off', 'BurstRejection','off', 'Distance','Euclidian' );
    
    % drop bad channels from EEG
   	rej_chans{ix} = setdiff({EEG.chanlocs.labels}, {ASR.chanlocs.labels});
    if numel(	rej_chans{ix} ) > 6
        warning('... >6 bad channels detected... please review dataset')
        continue
    end
    ASR.urchanlocs = EEG.chanlocs;

    % save corrected data
    ASR = pop_saveset( ASR, 'filename', [sname, '_asr.set'], 'filepath', proc_asr);

    % record number of channels rejected per subject
    fid = fopen( fname, 'a' );
    fprintf( fid, '%s, %d\n', sname, numel(rej_chans{ix}) );
    fclose( fid );

end

% save bad channel info
save(sprintf('WIM_HB_ASR_badChans_%s.mat', datestr(now, 'yyyy-mm-dd')), 'rej_chans')
