%%
clear
close

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))
addpath(genpath(path_heplab))
addpath(genpath(path_pncst))

%% batch process

HRVparams = WIM_HB_InitHRVparams;
manualScreen = 1;    % perform manual review of R-peak extraction (1 = yes, else skip)

subjs = dir(fullfile(wim_preproc, 'ecg', 'MWI*_ecg.set'));

for ix = 1:length(subjs) 
    
    sname = subjs(ix).name(1:6);
    fprintf(['Loading ' sname '...\n'])

    % load preprocessed ECG
    ECG = pop_loadset( 'filename', subjs(ix).name, 'filepath', fullfile(wim_preproc, 'ecg') );
    
    % z-normalise
    ECG.data = normalize(ECG.data);

    % annotate R-waves
    fprintf('...extracting QRS events...\n')
    [t,rr,jqrs_ann] = ConvertRawDataToRRIntervals( ECG.data, HRVparams, sname);
    ECG = updateEvents(ECG, num2cell(jqrs_ann));
    
    % review IBIs
%    checkIBI(ECG, 'QRS')

    % save annotated file
    ECG = pop_saveset( ECG, 'filename', [sname, '_ecg_qrs.set'], 'filepath', fullfile(wim_preproc, 'ecg') );

    
    % screen ECG for missing/incorrect R-peak annotations
    if manualScreen==1

        % setup HEP struct
        HEP.srate = ECG.srate;
        
        % get IBIs
        ibs = epochIBI(ECG);
        ibs2 = cell(size(ibs));

        for ex = 1:size(ibs,2)
    
            % heplab gui -- correct annotations & close window 
            % (do not press 'SAVE' -- function operates on EEGlab struct)
            HEP.sec_ini = 0;
            HEP.winsec = 20;
            HEP.ecg = heplab_ecg_filt(ECG.data(ibs{1,ex}(1)-ECG.srate/2:ibs{1,ex}(end)+ECG.srate/2), ECG.srate, 1, 20);
            HEP.qrs = ibs{1,ex}'-(ibs{1,ex}(1)-ECG.srate/2);
            heplab
    
            % await key press
            fprintf('Epoch %g/%g\n', ex, size(ibs,2))
            pause
    
            % record corrected peak latencies
            ibs2{1,ex} = HEP.qrs'+(ibs{1,ex}(1)-ECG.srate/2);
    
        end
            
        % update event struct with corrected events
        ECG.event = ECG.event(~strcmp({ECG.event.type}, 'QRS'));
        ECG = updateEvents(ECG, ibs2(1,:));

        % save corrected annotations
        ECG = pop_saveset( ECG, 'filename', [sname, '_ecg_qrs_screened.set'], 'filepath', fullfile(wim_preproc, 'ecg') );
    
    end

    
end

