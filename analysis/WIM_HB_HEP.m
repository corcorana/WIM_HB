%%
clear
close all

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))


%% batch process
subjs = dlmread(['..' filesep 'subjs.txt']);

% epoch/erp timing
ep_len = 30; % epoch duration (s);
cutoff = .25; % period between epoch offset & probe onset (s)
erp_on = -.5; % erp onset relative to R-peak trigger (s)
erp_off = 1; % erp offset relative to R-peak trigger (s)
min_ibi = 0; % minimal time to next R-peak (s)

[allERP, allRes] = deal(cell(length(subjs),1));

trl = [];
for ix = 1 :length(subjs)
    
    snum = subjs(ix);
    sname = ['MWI' num2str(snum)];
    fprintf(['Loading ' sname '...\n'])

    % load ICA-corrected EEG 
    EEG = pop_loadset( 'filename', [sname, '_ica_rej.set'], 'filepath', fullfile(eeg_preproc, 'ica'));
    fs = EEG.srate;

    % load corrected ECG data
    ECG = pop_loadset( 'filename', [sname, '_ibi_corrected_trimmed.set'], 'filepath', ['..' filesep 'ibi']);

    % load behav data 
    behav_name = dir([path_behav filesep '*s' num2str(snum) '*.mat']);
    load([path_behav filesep behav_name.name]);

    % import corrected event list and ECG data into EEG struct
    if EEG.times(end)~=ECG.times(end)
        warning('EEG and ECG record differ in length, skipping...')
        continue
    elseif sum(~strcmp({ECG.event.type}, 'ECG'))~=sum(~strcmp({ECG.event.type}, 'ECG'))
        warning('number of non-ECG event markers does not match, skipping...')
        continue
    else
        EEG.event = ECG.event;
        EEG.data = [EEG.data; ECG.data];
        EEG.nbchan = EEG.nbchan +1;
        [EEG.chanlocs(end+1).labels, EEG.chanlocs(end+1).type] = deal('ECG');
    end
    
    % low-pass
    EEG = pop_eegfiltnew(EEG, [], 30);

    % compute R-peak-locked ERPs in epochs preceding probes
    prbs = find(contains({EEG.event.type}, 'P  '))';
    eps = [[EEG.event(prbs).latency]'-ep_len*fs, [EEG.event(prbs).latency]'-cutoff*fs];
    rpks = [EEG.event(strcmp({EEG.event.type}, 'ECG')).latency]';

    t = erp_on:1/fs:erp_off;
    r = cell(length(prbs), 1);
    ERP = nan(EEG.nbchan, length(t), length(prbs));
    for px = 1:length(prbs)
        ridx = rpks > eps(px,1) & rpks < eps(px,2);
        ridx([false; diff(ridx)==-1])=true; % append final r peak to evaluate min_ibi
        r{px} = rpks(ridx);
        r{px}([diff(r{px}./fs) < min_ibi; true]) = nan; % exclude r peaks followed by short ibi
        if isempty(r{px})
            warning('no events in epoch %g', px)
            continue
        end
        trl = nan(EEG.nbchan, length(t), length(r{px}));
        for ex = find(~isnan(r{px})')
            trl(:,:,ex) = EEG.data(:,ceil(r{px}(ex)+fs*erp_on):ceil(r{px}(ex)+fs*erp_off));
        end
        ERP(:,:,px) = mean(trl,3, 'omitnan');
        ERP(:,:,px) = ERP(:,:,px) - median(ERP(:,t>=-.3&t<=-.1,px),2); 

    end

    % collate ERP & response data
    allERP{ix} = ERP;
    probe_res(probe_res(:,32)==4,32)=3;
    allRes{ix} = probe_res;
    
end

chanlocs = EEG.chanlocs;

save('WIM_HB_HEP.mat', 'allERP', 'allRes', 'ep_len', 'cutoff', 'min_ibi', 't', 'chanlocs')
