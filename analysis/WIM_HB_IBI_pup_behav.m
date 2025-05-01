%%
clear
close

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))
addpath(['..' filesep 'preproc'])

%% batch process
subjs = dir(fullfile(wim_preproc, 'ecg', '*_qrs_screened.set'));

% some analysis params
rfs = 5; % resample frequency (Hz) for timeseries interpolation
rt_thresh = 300; % RT cutoff (ms) -- exclude responses faster than this from resp speed estimates
ep_len = 10; % preprobe epoch length (s) -- period prior to probe onset

nexc = nan(length(subjs),2);
bn = nan(6,length(subjs));
tabs = [];
dtab = [];

for ix = 1 :length(subjs)
    
    snum = str2double(subjs(ix).name(4:6));
    sname = subjs(ix).name(1:6);
    fprintf(['Loading ' sname '...\n'])

    %% ECG data

    % load corrected ECG data
    ECG = pop_loadset( 'filename', subjs(ix).name, 'filepath', fullfile(wim_preproc, 'ecg') );
    fs = ECG.srate;

    % collect number of trials per block
    for bx = 1:6
        tmp = ECG.test_res(ECG.test_res(:,1)==bx,4);
        bn(bx,ix) = tmp(end);
    end

    % get RR interval latencies and calculate IBIs
    rri = epochIBI(ECG);
    ibi = cellfun( @(x){[nan, round(diff(x/fs)*1000)]}, rri);
    
    % reject improbable IBIs
    ibi2 = cell2mat(ibi);
    ibi2 = ibi2(ibi2>=300 & ibi2<=2000);  % exclude extreme IBIs from calculation
    muIBI = mean(ibi2); sdIBI = std(ibi2);
    id = cellfun( @(x){x >= max(300, muIBI-sdIBI*5) & x <= min(muIBI+sdIBI*5, 2000)}, ibi); % index valid IBIs
    ibs(1,:) = cellfun( @(x,idx) {[x(1) x(idx)]}, rri, id); % Latencies
    ibs(2,:) = cellfun( @(x,idx) {[x(1) x(idx)]}, ibi, id); % IBIs
    nexc(ix,:) = [sum(cell2mat(id)==0)-length(id), length(cell2mat(id))-length(id)];  % number of excluded IBIs / total IBIs

    % low-freq resampling for behav analysis
    tmz = cellfun(@(t){t/fs}, ibs(1,:)); 
    ibs(3,:) = cellfun(@(t,i) {interp1( t(2:end), i(2:end), t(2):1/rfs:t(end), 'spline')}, tmz(1,:), ibs(2,:) );  % ignore 1st timepoint/RRI (NaN)

    % compile summary statistics derived from preprobe epochs
    ep_off = [ECG.event(ismember({ECG.event.type}, {'S  3' 'P  1'})).latency];
    ep_on = ep_off-fs*ep_len;

    ii = nan(numel(ep_off),10);
    for ex = 1:size(ii,1)
        ii(ex,1:6) = [snum, ex, ECG.probe_res(ex, [4,5,32,38])]; % subj num, prb num, block num, stim type, mstate, vigil
        ii(ex,7) = mean( ibs{3,ex}(tmz{1,ex}(2):1/rfs:tmz{1,ex}(end)>ep_on(ex)/fs & tmz{1,ex}(2):1/rfs:tmz{1,ex}(end)<ep_off(ex)/fs) );
        ii(ex,8) = std( ibs{3,ex}(tmz{1,ex}(2):1/rfs:tmz{1,ex}(end)>ep_on(ex)/fs & tmz{1,ex}(2):1/rfs:tmz{1,ex}(end)<ep_off(ex)/fs)) / ...
            mean( ibs{3,ex}(tmz{1,ex}(2):1/rfs:tmz{1,ex}(end)>ep_on(ex)/fs & tmz{1,ex}(2):1/rfs:tmz{1,ex}(end)<ep_off(ex)/fs) );
        ii(ex,9) = mean( (ibs{3,ex}(tmz{1,ex}(2):1/rfs:tmz{1,ex}(end)>ep_on(ex)/fs & tmz{1,ex}(2):1/rfs:tmz{1,ex}(end)<ep_off(ex)/fs)-muIBI) / sdIBI );
        ii(ex,10) = sum(~isnan( ibs{2,ex}(ibs{1,ex}>ep_on(ex) & ibs{1,ex}<ep_off(ex)) ));
    end


    %% pupil data
    
    fname = fullfile(wim_preproc, 'eye', [sname, '_pup.mat']);

    if exist(fname, 'file')

        % load preprocessed pupil data
        load(fname)
        events = EL_events.Events;

        % epoched estimates
        pon = regexp(events.type, '^P[1-9]$|^P10$'); % events marking probe onset
        poff = regexp(events.type, '^EP[1-9]$|^EP10$'); % events marking probe offset
        tsta = events.time(~cellfun(@isempty, [{[]}, poff(1:end-1)])); % timelock to first marker following probe offset
        bsta = ~cellfun(@isempty, regexp(events.type, '^B1$')); % index block 1 onset
        tsta = [events.time([false, bsta(1:end-1)]), tsta(1:end-1)]; % include block 1 onset with restarts following probes
        tend = events.time(~cellfun(@isempty, pon)); % timelock to next probe onset
        
        pup = cell(1, length(tend));
        for ex = 1:length(pup) 
            tp = tsta(ex):tend(ex);
            pup{1,ex} = EL_data.time(EL_data.time>tsta(ex) & EL_data.time<tend(ex) )'; % latencies
            pup{2,ex} = EL_data.filt_pupilSize(EL_data.time>tsta(ex) & EL_data.time<tend(ex) )'; % filtered pup size
        end
        muPup = mean(cell2mat(pup(2,:)), 'omitnan'); sdPup = std(cell2mat(pup(2,:)), 'omitnan');
        pup(3,:) = cellfun( @(x) {(x-muPup)/sdPup}, pup(2,:) ); % z-normalise


        % compile summary statistics derived from preprobe epochs
        ep_on = tend-ep_len*EL_headers.Fs;
        ee = nan(numel(ep_on),3);
        for ex = 1:size(ee,1)
            ee(ex,1) = mean( pup{2,ex}(pup{1,ex}>=ep_on(ex)), 'omitnan' ); 
            ee(ex,2) = mean( pup{3,ex}(pup{1,ex}>=ep_on(ex)), 'omitnan' ); 
            ee(ex,3) = sum(~isnan(pup{3,ex}(pup{1,ex}>=ep_on(ex)))) / length(pup{3,ex}(pup{1,ex}>=ep_on(ex)));
        end

    else
        ee = nan(numel(ep_on),3);
    end
    
   
    %% behav data

    tlat = ECG.test_res(:,8); % trial times
    go = ECG.test_res(:, 12); % Go trials
    rt = ECG.test_res(:,10) - tlat; % response times
    rsp = rt;
    rsp(rsp<rt_thresh/1000 | isnan(go)) = nan;  % keep Go trial responses above RT threshold
    rsp = 1./rsp; % convert RTs to response speed

    % epoched estimates
    tend = num2cell(ECG.probe_res(:,2)); 
    tsta = num2cell([0; ECG.probe_res(1:end-1,3)]);

    rts(1,:) = cellfun(@(t1,t2) {tlat(tlat>t1 & tlat<t2)}, tsta, tend); % trl times
    rts(2,:) = cellfun(@(t1,t2) {rsp(tlat>t1 & tlat<t2) }, tsta, tend); % trimmed RS
    rts(3,:) = cellfun(@(t1,t2) {rt(tlat>t1 & tlat<t2) }, tsta, tend);  % uncorrected RT
    rts(4,:) = cellfun(@(p1,p2) {interp1(p1(~isnan(p2)), p2(~isnan(p2)), p1, 'linear')}, rts(1,:), rts(2,:) ); % interpolated RS
    rts(5,:) = cellfun(@(p1,p2) {interp1(p1(~isnan(p2)), p2(~isnan(p2)), p1, 'linear')}, rts(1,:), rts(3,:) ); % interpolated RT

    % mark non-responses on NoGo invalid if neighboured by non-responses on Go
    ng = ECG.test_res(:, 11); % NoGo trials
    ngx = [ng, [0; abs(go(1:end-1)-1)], abs([go(2:end)-1; 0]) ]; % compare with lag / lead Go trial response
    ngx([false; diff(tlat)>1.3],2) = nan; % exclude preceding trials from previous epoch
    ngx([diff(tlat)>1.3;false],3) = nan; % exclude subsequent trials from next epoch
    ng(sum(ngx,2)==3) = nan;
    
    % compile summary statistics derived from preprobe epochs
    ep_off = ECG.probe_res(:,2);
    ep_on = ep_off-ep_len;
    
    tt = nan(numel(ep_off),7);
    for ex = 1:size(tt,1)
        % response speed
        tt(ex,1) = mean( rts{2,ex}(rts{1,ex}>ep_on(ex) & rts{1,ex}<ep_off(ex)), 'omitnan' ); 
        tt(ex,2) = (std( rts{4,ex}(rts{1,ex}>ep_on(ex) & rts{1,ex}<ep_off(ex)), 'omitnan' ) / tt(ex,1) );
        tt(ex,3) = sum(~isnan( rts{2,ex}(rts{1,ex}>ep_on(ex) & rts{1,ex}<ep_off(ex)) ));
        % response accuracy
        tt(ex,4) = sum(ng(tlat>ep_on(ex) & tlat<ep_off(ex))==1); % nogo correct ("CR")
        tt(ex,5) = sum(ng(tlat>ep_on(ex) & tlat<ep_off(ex))==0); % nogo error ("FA") 
        tt(ex,6) = sum(go(tlat>ep_on(ex) & tlat<ep_off(ex))==1); % go correct ("H")
        tt(ex,7) = sum(go(tlat>ep_on(ex) & tlat<ep_off(ex))==0); % go error ("M") 
    end

    % collate tables
    tabs = [tabs; ii, ee, tt];


    %% cardiac deceleration

    prb_on = [ECG.event(ismember({ECG.event.type}, {'S  3' 'P  1'})).latency];

    dd = [];
    for ex = 1:length(prb_on)
        
        % isolate preprobe NoGo trial onset times and response data
        tep = ECG.test_res(tlat>=ep_on(ex)&tlat<ep_off(ex), 8);
        rep = ECG.test_res(tlat>=ep_on(ex)&tlat<ep_off(ex), 11);
        tep = tep(~isnan(rep));
        rep = rep(~isnan(rep));
    
        % allow for >1 NoGo trial per epoch
        for nx = 1:numel(rep)
            % collate subject, probe, stim, state, vigil, response, time
            pinfo = [snum, ex, ECG.probe_res(ex, [5, 32,38]), rep(nx), tep(nx)]; 

            % localise IBI containing NoGo trial onset
            tims = ep_off(ex) - tep(nx); % time interval between NoGo trial and probe onset
            ib0 = prb_on(ex) - tims*fs;  % latency between NoGo trial and probe onset
            rpk = find(ibs{1,ex}>ib0); % R peaks intervening between NoGo trial and probe onset

            % obtain IBI series
            try
                d = ibs{2,ex}(rpk(1)-3:rpk(1)+3);
            catch
                try
                    d = [ibs{2,ex}(rpk(1)-3:rpk(1)+2), nan];
                catch
                    try
                        d = [ibs{2,ex}(rpk(1)-3:rpk(1)+1), nan, nan];
                    catch
                        % abort if no IBIs post IBI-0
                    end
                end
            end
            dd = [dd; [pinfo, d] ];
        end
    end
    dd(:,8:14) = (dd(:,8:14)-muIBI)./sdIBI; % z-normalise IBIs
    
    % collate across subjects
    dtab = [dtab; dd];

end


%% output tables

% epoch means
varlabs = {'subj_id', 'probe_num', 'block_num', 'stim_type', 'state', 'vigil', ...
    'muIBI', 'cvIBI', 'zuIBI', 'nIBI', 'muPup', 'zuPup', 'prPup', 'muRS', 'cvRS', 'nRS', 'CR', 'FA', 'H', 'M' };
tab = array2table(tabs, 'VariableNames', varlabs);

% recover 3 digit subject ID code
tab.subj_id = num2str(tab.subj_id, '%03.f');

% recode stim type
tab.stim_type(tab.stim_type==1)='F';
tab.stim_type(tab.stim_type==2)='D';
tab.stim_type = char(tab.stim_type);

% output
writetable(tab, 'WIM_HB_IBI_pup_behav.csv')


% cardiac deceleration
varlabs = {'subj_id', 'probe_num', 'stim_type', 'state', 'vigil', 'resp_type', 'trl_time', ...
    'IBI_3', 'IBI_2', 'IBI_1', 'IBI_0', 'IBI.1', 'IBI.2', 'IBI.3' };
tab = array2table(dtab, 'VariableNames', varlabs);

tab.stim_type(tab.stim_type==1)='F';
tab.stim_type(tab.stim_type==2)='D';
tab.stim_type = char(tab.stim_type);

tab.resp_type(tab.resp_type==1)='C';
tab.resp_type(tab.resp_type==0)='F';
tab.resp_type = char(tab.resp_type);

% output
writetable(tab, 'WIM_HB_IBI_NoGo.csv' )
