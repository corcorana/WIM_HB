%%
clear
close

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))

proc_ica = fullfile(wim_preproc, 'ica');

%% batch process

subjs = dir(fullfile(wim_preproc, 'ica', 'MWI*_ica_rej.set'));

fname = sprintf('WIM_HB_ICA_rejCFA_%s.csv', datestr(now, 'yyyy-mm-dd'));
fid = fopen( fname, 'w' );
fprintf( fid, '%s, %s\n', 'Subject', 'rejCFA');
fclose( fid );

for ix = 1:length(subjs)

    sname = subjs(ix).name(1:6);
    fprintf(['Loading ' sname '...\n'])

    % load ICA'd data
    EEG = pop_loadset('filename', [sname, '_ica_rej.set'], 'filepath', proc_ica );
    
    % epoch around R-peak
    CFA = pop_epoch( EEG, { 'QRS' }, [-0.2 0.2], 'epochinfo', 'yes'); 

    % run ICA
    CFA =  pop_runica(CFA, 'icatype', 'runica', 'extended',1, 'interrupt','off');
    
    % calculate component-wise mean coherence 
    phidiff = nan(CFA.pnts,CFA.trials);
    cohere = nan(size(CFA.icaact,1),1);
    lpfrq = 25;
    for cx = 1: size(CFA.icaact, 1)
        fprintf('... component #%g/%g\n', cx, size(CFA.icaact, 1))
        for tx = 1:CFA.trials
            x = lowpass(CFA.icaact(cx,:,tx), CFA.srate, lpfrq);
            ecg = lowpass(CFA.data(strcmp({CFA.chanlocs.type}, 'ECG'),:,tx), CFA.srate, lpfrq);
            phidiff(:,tx) = angle(hilbert(x)) - angle(hilbert(ecg));
        end
        cohere(cx) = mean( abs(mean(exp(1i*phidiff),2)) );
    end

    % identify up to nmax high-coherence IC-ECG pairs
    nmax = 3;
    comprej = nan(1,nmax);
    for rx = 1:length(comprej)
        tmp = find(cohere>(mean(cohere, 'omitnan')+std(cohere, 'omitnan')*3));
        if isempty(tmp)
            break
        else
            [B, I] = sort(cohere(tmp), 'descend');
            comprej(rx) = tmp(I(1));
            cohere(tmp(I(1))) = nan;
        end
    end

    % subtract from original dataset
    EEG.icaact = [];
    EEG.icawinv = CFA.icawinv;
    EEG.icasphere = CFA.icasphere;
    EEG.icaweights = CFA.icaweights;
    
    comprej = comprej(~isnan(comprej));
    EEG = pop_subcomp(EEG, comprej);

    % record number of rejected components per subject
    fid = fopen( fname, 'a' );
    fprintf( fid, '%s, %d\n', sname, numel(comprej) );
    fclose( fid );

    % save decomposition & corrected data
    EEG = pop_saveset( EEG, 'filename', [sname, '_ica_rej_cfa.set'], 'filepath', proc_ica );

end



%% lowpass filter
function LP = lowpass(timecourse, SamplingRate, f_cut, filterOrder)
% custom function adapted from https://github.com/andrillon/LSCPtools/maths

if (nargin < 4)
    filterOrder = 2;
end

[b, a] = butter(filterOrder, (f_cut/SamplingRate)*2,'low');
LP = filtfilt(b, a, timecourse );

 
end
