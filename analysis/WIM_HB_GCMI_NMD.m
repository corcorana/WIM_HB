%%
clear
close all

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))
addpath(genpath(path_nmd))
addpath(genpath(path_gcmi))


%% params

dfs = 50;   % downsampled frequency (Hz)
eplen = 10; % probe epoch analysis window length (s)
wlen = 16;  % window length (segment prior to probe; s)

% timelags (s) 
hbLags = logspace(0, 1.1, 12)/10;   % heart to brain 
bhLags = logspace(0.2, 1.4, 12)/10; % brain to heart (slower response at SA node)

% frequency bands (Hz)
cfrex = logspace(0, 1.2, 8)'; 
fbands = [cfrex*.8, cfrex*1.2];

% gcmi bias correction
bc = false;


%% batch process
subjs = dir(fullfile(wim_preproc, 'bhi', 'MWI*_bhi.set'));

MI = cell(length(subjs),2);
for ix = 1 :length(subjs)
    
    sname = subjs(ix).name(1:6);
    fprintf(['Loading ' sname '...\n'])
   
    % load processed EEG
    EEG = pop_loadset( 'filename', [sname, '_bhi.set'], 'filepath', fullfile(wim_preproc, 'bhi'));
    
    % downsample
    EEG = pop_resample( EEG, dfs );

    % separate ECG
    ECG = pop_select( EEG, 'channel', {'ECG'});
    EEG = pop_select( EEG, 'rmchantype', {'ECG'});

    % epoch ECG
    ECG = pop_epoch( ECG, { 'S  3', 'P  1' }, [-wlen, 2], 'epochinfo', 'yes');

    % compute copnorm'd ECG phi for each epoch
    cibi_phi = nan(ECG.pnts-2*dfs,2,ECG.trials);
    for nEp = 1:ECG.trials
        [tfr, freq, wopt] = wt( ECG.data(1,:,nEp), dfs, 'fmin', .5, 'fmax', 2, 'Plot', 'off', 'Display', 'off');
        tfSupp = ecurve(tfr, freq, wopt);
        [~, ~, ~, asig] = rectfr_mod(tfSupp(1,:), tfr, freq, wopt); % function adapted to output 'asig'
        ibi_phi = [ real(asig)./abs(asig); imag(asig)./abs(asig) ]';
        cn = copnorm(ibi_phi); 
        cibi_phi(:,:,nEp) = cn(1:end-2*dfs,:);
    end

    % EEG frequency band
    for nFrq = 1 : size(fbands, 1)
     
        % filter freq band 
        EEGband = pop_eegfiltnew(EEG, fbands(nFrq, 1), []);
        EEGband = pop_eegfiltnew(EEGband, [], fbands(nFrq, 2));

        % epoch
        EEGband = pop_epoch( EEGband, { 'S  3', 'P  1' }, [-wlen, 2], 'epochinfo', 'yes');

        for nEp = 1:EEGband.trials

            % compute copnorm'd EEG phi for each channel & epoch
            hilb_eeg = hilbert( squeeze(EEGband.data(:,:,nEp)') );  
            eeg_amp = abs( hilb_eeg );
            eeg_phi = permute( cat( 3, real(hilb_eeg)./abs(hilb_eeg), imag(hilb_eeg)./abs(hilb_eeg) ), [1,3,2] );
            cn = copnorm(eeg_phi);
            cn(isnan(eeg_phi)) = nan;
            ceeg_phi = cn(1:end-2*dfs,:,:);

            % heart -> brain MI (lead ibi)
            deeg = ceeg_phi(((wlen-eplen)*dfs)+1:end,:,:); % fix at lag-0
            for nLag = 1:length(hbLags)
                d = round(hbLags(nLag)*dfs);
                dibi = cibi_phi((wlen-eplen)*dfs-d+1:(end-d),:,nEp);
                for nCh = 1:EEG.nbchan
                    try
                        MI{ix,1}(1,nCh,nFrq,nLag,nEp) = mi_gg(deeg(:,:,nCh), dibi, bc, true);
                    catch
                        MI{ix,1}(1,nCh,nFrq,nLag,nEp) = nan;
                    end
                end
            end
            
            % brain -> heart MI (lead eeg)                            
            dibi = cibi_phi(((wlen-eplen)*dfs)+1:end,:,nEp); % fix at lag-0
            for nLag = 1:length(bhLags)
                d = round(bhLags(nLag)*dfs);
                deeg = ceeg_phi((wlen-eplen)*dfs-d+1:(end-d),:,:);
                for nCh = 1:EEG.nbchan
                    try
                        MI{ix,1}(2,nCh,nFrq,nLag,nEp) = mi_gg(dibi, deeg(:,:,nCh), bc, true);
                    catch
                        MI{ix,1}(2,nCh,nFrq,nLag,nEp) = nan;
                    end
                end
            end

        end
    end

    % collate probe info
    EEG.probe_res(EEG.probe_res(:,32)>3,32)=3;
    MI{ix,2} = [ repmat(str2double(sname(4:6)), size(EEG.probe_res,1), 1), [1:size(EEG.probe_res,1)]', EEG.probe_res(:, [5,32,38])];

end

chanlocs = EEG.chanlocs;

save('WIM_HB_GCMI_NMD.mat', 'MI', 'eplen', 'dfs', 'hbLags', 'bhLags', 'fbands', 'cfrex', 'chanlocs')
