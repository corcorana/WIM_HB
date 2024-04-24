%%
clear
close all

%% paths
run(['..' filesep 'localdef_WIM_ECG'])
addpath(genpath(path_eeglab))
addpath(genpath(path_heplab))

proc_ibi = fullfile('..', filesep, 'ibi');
if ~exist(proc_ibi, 'dir')
    mkdir(proc_ibi)
end

%% batch process
subjs = dlmread(['..' filesep 'subjs.txt']);

for ix = 1:length(subjs)
    
    snum = subjs(ix);
    sname = ['MWI' num2str(snum)];
    fprintf(['Loading ' sname '...\n'])

    % load trimmed & annotated ECG data
    EEG = pop_loadset( 'filename', [sname, '_min.set'], 'filepath', fullfile(eeg_preproc, 'min'));
    ECG = pop_select( EEG, 'channel', {'ECG'});
    fs = ECG.srate;
    
    % setup HEP struct
    HEP.srate = fs;
    
    % get IBI latencies
    ibs = indexIBI(ECG);
    ibs2 = cell(size(ibs));

    % screen ECG for missing/incorrect R-peak annotations
    for ex = 1:size(ibs,2)

        % heplab gui -- correct annotations & close window 
        % (do not press 'SAVE' -- function operates on EEGlab struct)
        HEP.sec_ini = 0;
        HEP.winsec = 20;
        HEP.ecg = heplab_ecg_filt(ECG.data(ibs{1,ex}(1)-fs/2:ibs{1,ex}(end)+fs/2), ECG.srate, 1, 20);
        HEP.qrs = ibs{1,ex}'-(ibs{1,ex}(1)-fs/2);
        heplab

        % await key press
        fprintf('Epoch %g/%g\n', ex, size(ibs,2))
        pause

        % record corrected peak latencies
        ibs2{1,ex} = HEP.qrs'+(ibs{1,ex}(1)-fs/2);

    end
        
    % update event struct & review
    ECG2 = correctRRI(ECG, ibs2);
    checkRRI(ibs2, ECG2)

    % save corrected file
    ECG2 = pop_saveset( ECG2, 'filename', [sname, '_ibi_corrected.set'], 'filepath', proc_ibi );
    

end
       