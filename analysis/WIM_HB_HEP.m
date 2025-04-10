%%
clear
close

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))


%% batch process
subjs = dir(fullfile(wim_preproc, 'bhi', 'MWI*_bhi.set'));

% epoch/erp timing
ep_len = 10;  % epoch duration (s);
cutoff = .3;  % period between epoch offset & probe onset (s)
erp_on = -.3; % erp onset relative to R-peak trigger (s)
erp_off = .9; % erp offset relative to R-peak trigger (s)
min_ibi = .6; % minimum time to next R-peak (s)
min_rpk = 5;  % minimum number of R-peaks in epoch to estimate preprobe HEP
lpf_cut = 20; % low-pass filter cutoff (Hz)

allHEP = cell(length(subjs), 2);

trl = [];
for ix = 1 :length(subjs)
    
    sname = subjs(ix).name(1:6);
    fprintf(['Loading ' sname '...\n'])

    % load processed EEG 
    EEG = pop_loadset( 'filename', [sname, '_bhi.set'], 'filepath', fullfile(wim_preproc, 'bhi'));

    % low-pass filter
    EEG = pop_eegfiltnew(EEG, [], lpf_cut);
    dat = EEG.data;
    fs = EEG.srate;

    % compute R-peak-locked ERPs in epochs preceding probes
    prbs = find(ismember({EEG.event.type}, {'S  3' 'P  1'})); % probe onsets
    eps = [[EEG.event(prbs).latency]'-ep_len*fs, [EEG.event(prbs).latency]'-cutoff*fs]; % epoch bounds
    rpks = [EEG.event(strcmp({EEG.event.type}, 'QRS')).latency]';

    t = erp_on:1/fs:erp_off;
    HEP = nan(EEG.nbchan, length(t), length(prbs));
    for px = 1:length(prbs)
        ridx = rpks > eps(px,1) & rpks < eps(px,2);
        ridx([false; diff(ridx)==-1])=true; % append final R peak to evaluate min_ibi
        r = rpks(ridx);
        if isempty(r)
            warning('no events in epoch %g', px)
            continue
        end
        r([diff(r./fs) < min_ibi; true]) = nan; % exclude R peaks followed by short ibi (& appended R peak)
        if sum(~isnan(r)) < min_rpk
            warning('insufficent R-peaks, skipping epoch %g', px)
            continue
        else
            trl = nan(EEG.nbchan, length(t), sum(~isnan(r)));
            for rx = find(~isnan(r)')
                trl(:,:,rx) = dat(:,ceil(r(rx)+fs*erp_on):ceil(r(rx)+fs*erp_off));
            end
            HEP(:,:,px) = mean(trl,3, 'omitnan');
        end
    end

    % collate ERP & response data
    allHEP{ix,1} = HEP;
    EEG.probe_res(EEG.probe_res(:,32)>3,32)=3;
    allHEP{ix,2} = [ repmat(str2double(sname(4:6)), size(EEG.probe_res,1), 1), [1:size(EEG.probe_res,1)]', EEG.probe_res(:, [5,32,38])];

end

chanlocs = EEG.chanlocs;

save('WIM_HB_HEP.mat', 'allHEP', 'ep_len', 'erp_on', 'erp_off', 'cutoff', 'min_ibi', 'min_rpk', 'lpf_cut', 't', 'fs', 'chanlocs')
