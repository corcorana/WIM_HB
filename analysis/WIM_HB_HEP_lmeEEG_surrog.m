%%
clear
close

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))
addpath(genpath(path_lmeEEG))

%% batch process
subjs = dir(fullfile(wim_preproc, 'bhi', 'MWI*_bhi.set'));

load('WIM_HB_HEP.mat')
load('WIM_HB_HEP_lmeEEG.mat')

%% import & process data to be permuted

DAT = cell(length(subjs),2);
for ix = 1:length(subjs)
           
    sname = subjs(ix).name(1:6);
    fprintf(['Loading ' sname '...\n'])
     
    % load processed EEG 
    EEG = pop_loadset( 'filename', [sname, '_bhi.set'], 'filepath', fullfile(wim_preproc, 'bhi'));
        
    % low-pass filter
    EEG = pop_eegfiltnew(EEG, [], lpf_cut);
    
    % add to array
    DAT(ix,:) = {EEG.data, EEG.event};

end


%% permutations (MS only)
seg_len = 40; % length of segment prior to probe for R-peak permutation (sec)
nperms = 100; % number of times R-peak series permuted
t_surr = nan(nperms, length(chanlocs), length(t), length(mod.CoefficientNames) ); % Initialize 

for np = 1:nperms

    allSUR = cell(length(subjs), 2);

    for ix = 1:size(DAT,1)

        % compute R-peak-locked ERPs in epochs preceding probes
        eeg = DAT{ix,1};
        evt = DAT{ix,2};
        prbs = find(ismember({evt.type}, {'S  3' 'P  1'})); % probe onsets
        rpks = [evt(strcmp({evt.type}, 'QRS')).latency];

        sHEP = nan(size(eeg,1), length(t), length(prbs));
        for px = 1:length(prbs)  
            seg_on = evt(prbs(px)).latency-(seg_len*fs);
            seg_off = evt(prbs(px)).latency;
            ridx = rpks > seg_on & rpks < seg_off;
            ibs = diff(rpks(ridx));
            rp = randperm(length(ibs));
            spks = cumsum([seg_on + round(rand(1)*fs), ibs(rp)]);
            sidx = spks > seg_off-(ep_len*fs) & spks < seg_off-(cutoff*fs);
            s = spks(sidx);
            s(diff(s./fs) < min_ibi) = nan; % exclude r peaks followed by short ibi
            if sum(~isnan(s)) < min_rpk
                warning('insufficent R-peaks, skipping epoch %g', px)
                continue
            else
                trl = nan(size(eeg,1), length(t), length(s));
                for sx = find(~isnan(s))
                    trl = eeg(:,ceil(s(sx)+fs*erp_on):ceil(s(sx)+fs*erp_off));
                end
                sHEP(:,:,px) = mean(trl,3, 'omitnan');
            end
        end
            
        % append to surrogate HEP matrix
        allSUR{ix,1} = sHEP;

    end

    sDAT = cat(3, allSUR{:,1});
        
    % define state variable & exclude epochs with missing data
    u = unique(TAB.sid(TAB.ms~="ON")); % identify subjects that report OFF states
    eps = ismember(TAB.sid, u) & ~isundefined(TAB.ms); % exclude subjects that don't report OFF; any missing probe reports
    x = TAB.ms(eps);
    tab = [TAB(eps,:), table(x)];
    dat = sDAT(:,:,eps);

    % fit mixed models & perform mass univariate regression
    fprintf('Obtaining surrogate statistics [%g/%g]...\n', np, nperms)
    t_surr(np,:,:,:) = lmeEEG_fitMods(dat, tab, mform);

end

Results(1).t_surr = t_surr;

% save results
save('WIM_HB_HEP_lmeEEG.mat', 'Results', '-append')
