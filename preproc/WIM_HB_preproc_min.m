%%
clear all
close all

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))
addpath(genpath(path_pncst))

raw_eeg = dir([path_raweeg filesep '*.vhdr']);

proc_min = fullfile(eeg_preproc, 'min');
if ~exist(proc_min, 'dir')
    mkdir(proc_min)
end

%% batch process

HRVparams = WIM_HB_InitHRVparams;

subjs = dlmread(['..' filesep 'subjs.txt']);

for ix = 1:length(subjs)
    snum = subjs(ix);
    sname = ['MWI' num2str(snum)];
    fname = raw_eeg(contains({raw_eeg.name}, sname)).name;
    fprintf(['Loading ' sname '...\n'])

    behav_name = dir([path_behav filesep '*s' num2str(snum) '*.mat']);
    if isempty(behav_name)
        fprintf('... Sub %g SKIPPED: no behav file\n', snum)
        fid = fopen('warnings.txt','a');
        fprintf(fid,'%s ', sname, ': unable to locate behav .mat file');
        fprintf(fid,'\n');
        fclose(fid);
        continue
    else
        load([path_behav filesep behav_name.name]);
    end
    
    if size(probe_res, 1) ~= 60
        warning('Unexpected number of probe responses, skipping...\n')
        fid = fopen('warnings.txt','a');
        fprintf(fid,'%s ', sname, ': num probe res ==', num2str(size(probe_res, 1)));
        fprintf(fid,'\n');
        fclose(fid);
        continue
    end
    
    % load raw EEG
    EEG = pop_fileio( fullfile( path_raweeg, fname ), 'dataformat', 'brainvision_eeg' );

    % check sample rate
    fs = EEG.srate;
    if fs ~= 500
        warning('Unexpected sample rate (%g Hz), skipping...\n', fs)
        fid = fopen('warnings.txt','a');
        fprintf(fid,'%s ', sname, ': srate ==', num2str(fs) );
        fprintf(fid,'\n');
        fclose(fid);
        continue
    end
    
    % drop unused channels
    EEG = pop_select( EEG, 'nochannel',{'EOG' 'D1'});    
    if EEG.nbchan ~= 64
        warning( 'Unexpected number of EEG channels (n = %g), skipping...\n', size(EEG.data, 1) )
        fid = fopen('warnings.txt','a');
        fprintf(fid,'%s ', sname, ': num EEG chans ==', num2str(EEG.nbchan) );
        fprintf(fid,'\n');
        fclose(fid);
        continue
    end
    
    % trim off edges
    esta = find(ismember({EEG.event.type}, 'B  1'));
    esta = esta(end-5);
    eend = find(ismember({EEG.event.type}, 'K  1'));
    eend = eend(end);
    EEG = pop_select( EEG, 'point',[EEG.event(esta).latency, EEG.event(eend).latency+fs*5] );
    
    if length(find(ismember({EEG.event.type}, 'B  1'))) ~= 6
        warning('Unexpected number of blocks, skipping...\n')
        fid = fopen('warnings.txt','a');
        fprintf(fid,'%s ', sname, ': num blocks ==', num2str(length(find(ismember({EEG.event.type}, 'B  1')))));
        fprintf(fid,'\n');
        fclose(fid);
        continue
    end
    
    % index probes
    pidx = find(ismember({EEG.event.type}, 'P  1'));
    if length(pidx) ~= 60
        warning('Unexpected number of probe markers, skipping...\n')
        fid = fopen('warnings.txt','a');
        fprintf(fid,'%s ', sname, ': num probe markers ==', num2str(length(pidx)) );
        fprintf(fid,'\n');
        fclose(fid);
        continue
    end 
    
    % recode probe triggers according to reported mental state
    for px = 1:length(pidx)
        mx = probe_res(px,32);
        EEG = pop_editeventvals(EEG, 'changefield', {pidx(px), 'type', sprintf('P  %g%g', mx, px)} );
    end
    
    % index block trials -- assumes last trial marker is not a trial onset
    bsta = find(ismember({EEG.event.type}, 'B  1'));
    bend = find(ismember({EEG.event.type}, 'K  1'));
    tidx = find(ismember({EEG.event.type}, 'T  1'));    
    for bx = 1:length(bsta)
        fprintf('Finding trial markers, block %g\n', bx)
        trl = tidx(tidx > bsta(bx) & tidx < bend(bx));
        block_res = test_res(test_res(:,1)==bx,:);
        if block_res(1, 2) == 1     % face stimuli
            sx = 'F';
        elseif block_res(1, 2) == 2 % digit stimuli
            sx = 'D';
        end
        if length(trl)-1 == size(block_res,1)
            for nx = 1:length(trl)-1
                if ~isnan(block_res(nx, 12))
                    cx = 'G';
                elseif ~isnan(block_res(nx, 11))
                    cx = 'N';
                else
                    warning('Unable to resolve condition, trial %g\n', nx)
                    fid = fopen('warnings.txt','a');
                    fprintf(fid,'%s ', sname, ': problem resolving trial condition, block %g trial %g', bx, nx);
                    fprintf(fid,'\n');
                    fclose(fid);
                end
                EEG = pop_editeventvals(EEG, 'changefield', {trl(nx), 'type',...
                    sprintf('T  %g%s%s%g', bx, sx, cx, nx)} );
            end
        else
            warning('Unable to match trial markers and numbers, skipping...\n')
            fid = fopen('warnings.txt','a');
            fprintf(fid,'%s ', sname, ': problem matching trials from block', num2str(bx));
            fprintf(fid,'\n');
            fclose(fid);
            continue
        end
    end
    % insert response markers (all trials recoded as events are resorted by call to
    % pop_editeventvals above, invalidating order info for bx>1)
    for bx = 1:length(bsta)
        fprintf('Finding response markers, block %g\n', bx)
        trl = tidx(tidx > bsta(bx) & tidx < bend(bx));
        block_res = test_res(test_res(:,1)==bx,:);
        k = size(EEG.event,2); % number of existing events in the EEG structure
        urevents = num2cell( max([EEG.event.urevent])+1:max([EEG.event.urevent])+size(block_res,1) );
        resps = [EEG.event(trl(1:end-1)).latency]' + (block_res(:,10)-block_res(:,8))*fs;
        code = block_res(~isnan(resps), 11)==0;
        evt = num2cell( resps(~isnan(resps)) );
        types = repmat({'R  1'},1, length(evt));
        types(code) = deal({'R  0'});   % mark failure to inhibit response on NoGo
        durations = num2cell(ones(1, length(evt)));
        [EEG.event(1,k+1:k+length(evt)).latency] = evt{:}; % assign latencies
        [EEG.event(1,k+1:k+length(evt)).type] = types{:}; % assign types
        [EEG.event(1,k+1:k+length(evt)).duration] = durations{:}; % assign durations
        [EEG.event(1,k+1:k+length(evt)).urevent] = urevents{:}; % assign event index
    end
    
    % import channel locations and specify online ref (AFz)
    EEG = pop_chanedit(EEG, 'append', EEG.nbchan, 'changefield',{EEG.nbchan+1, 'labels' 'AFz'},...
        'lookup', fullfile( path_dipfit.folder, path_dipfit.name, 'standard_BESA','standard-10-5-cap385.elp'),...
        'setref', {find(~contains({EEG.chanlocs.labels}, 'dir')), 'AFz'});

    % filter copy of ECG for R-wave identification
    ECG = pop_select( EEG, 'channel',{'ECG'});     
    ECG = pop_eegfiltnew(ECG, 1, []);
    ECG = pop_eegfiltnew(ECG, [], 40);
    
    % identify R-wave peaks
    fprintf('...extracting QRS events...\n')
    [t,rr,jqrs_ann,SQIjw,StartSQIwindows_jw] = ConvertRawDataToRRIntervals(ECG.data/1000, HRVparams, sname);
    qrs = jqrs_ann;

    % transfer events to EEG.event (ensure unique urevents)
    k = size(EEG.event,2);
    urevents = num2cell( max([EEG.event.urevent])+1:max([EEG.event.urevent])+length(qrs) );
    evt = num2cell(qrs);
    types = repmat({'ECG'},1, length(evt));
    durations = num2cell(ones(1, length(evt)));
    [EEG.event(1,k+1:k+length(evt)).latency] = evt{:}; % assign latencies
    [EEG.event(1,k+1:k+length(evt)).type] = types{:}; % assign types
    [EEG.event(1,k+1:k+length(evt)).duration] = durations{:}; % assign durations
    [EEG.event(1,k+1:k+length(evt)).urevent] = urevents{:}; % assign event index

    % sort events by latency
    [~, idx] = sort([EEG.event.latency]);
    EEG.event = EEG.event(idx);

    % save .set files
    EEG = pop_saveset( EEG, 'filename', [sname, '_min.set'], 'filepath', proc_min);

end
