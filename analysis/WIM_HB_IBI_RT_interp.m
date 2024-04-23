%%
clear
close all

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))
addpath( ['..' filesep 'preproc'] )

%% batch process
subjs = dlmread(['..' filesep 'subjs.txt']);

% condition labels
msCat = {'ON', 'MW', 'MB'};

% some analysis params
rs = 5; % resample frequency (Hz) for timeseries interpolation
dfs = 50; % resample frequency (Hz) for interpolation to use in HBI analysis
rt_thresh = 300; % RT cutoff (ms) -- exclude responses faster than this from resp speed estimates
len = 30*rs; % pre-probe epoch duration (sample points)


ibi_tab = [];
trl_tab = [];
interpt = [];
for ix = 1 :length(subjs)
    
    snum = subjs(ix);
    sname = ['MWI' num2str(snum)];
    fprintf(['Loading ' sname '...\n'])

    % load corrected ECG data
    ECG = pop_loadset( 'filename', [sname, '_ibi_corrected.set'], 'filepath', ['..' filesep 'ibi'] );
    fs = ECG.srate;
    
    % load behav data (for probes)
    behav_name = dir([path_behav filesep '*s' num2str(snum) '*.mat']);
    load([path_behav filesep behav_name.name]);

    % get RR interval latencies and calculate IBIs
    rri = indexIBI(ECG);
    ibi = cellfun( @(x){[nan, round(diff(x/fs)*1000)]}, rri);
    
    % reject improbable IBIs
    ibi2 = cell2mat(ibi);
    ibi2 = ibi2(ibi2>=300&ibi2<=2000);  % exclude extreme IBIs
    muIBI = mean(ibi2); sdIBI = std(ibi2);
    id = cellfun( @(x){x >= muIBI-sdIBI*5 & x <= muIBI+sdIBI*5}, ibi); % index valid IBIs
    ibs(1,:) = cellfun( @(x,idx) {[x(1) x(idx)]}, rri, id); % Latencies
    ibs(2,:) = cellfun( @(x,idx) {[x(1) x(idx)]}, ibi, id); % IBIs

    % demean & normalise IBIs to SD of all (cleaned) epochs prior to interpolation
    tmz = cellfun(@(t){t/fs}, ibs(1,:)); 

    % low-freq resampling for behav analysis
    ibs(3,:) = cellfun(@(t,i) {interp1( t(2:end), i(2:end), t(2):1/rs:t(end), 'spline')}, tmz(1,:), ibs(2,:) );  % ignore 1st timepoint/RRI (NaN)
    
    % high-freq resampling for HBI analysis
    ibs(4,:) = cellfun(@(t,i) {interp1( t(2:end), i(2:end), t(2):1/dfs:t(end), 'spline')}, tmz(1,:), ibs(2,:) );  % ignore 1st timepoint/RRI (NaN)

    % get probe latencies
    prbs = {ECG.event(contains({ECG.event.type}, 'P')).latency}; 

    % collate ibi table
    ii = [];
    for px = 1:length(prbs)  
        ii = [ii; [ repmat([snum, px, probe_res(px,[1,4,5,32,38])], numel(ibs{1,px}), 1), ...
            ibs{1,px}'/fs, (ibs{1,px}-prbs{px})'/fs, reshape([ibs{2,px}], numel(ibs{1,px}), []) ] ];
    end
    ibi_tab = [ ibi_tab; [ii, normalize(ii(:,10))] ];


    %% RT data from trial
   
    % extract latencies
    rlat = [ECG.event(contains({ECG.event.type}, {'R'})).latency];
    t = regexp({ECG.event.type},'T  (\d)([FD])', 'tokens');
    tlat = [ECG.event(~cellfun(@isempty, t)).latency];
    
    % for all response markers, locate nearest preceding trial onset
    rmat = round((rlat-tlat')/fs*1000);     % col = responses, row = trial onsets
    rmat(rmat<0) = nan;   
    [rt, rti] = min(rmat); 
    
    % convert RTs to response speed (excluding fast responses)
    rsp = rt;
    rsp(rsp<rt_thresh)=nan;
    rsp = 1000./rsp;

    
    % impute resps amongst no-resp trials
    nr = nan(2, length(tlat));
    nr(1,rti) = rsp;
    nr(2,rti) = rt;
    
    % exclude NoGo responses 
    go_resp = test_res(:, 12); % index responses on Go trials
    nr(:,isnan(go_resp)) = nan; % nan any estimates on NoGo responses

    prbs1 = [{0}, prbs(1:end-1)];
    
    rts(1,:) = cellfun(@(p1,p2) {tlat(tlat>p1 & tlat<=p2)}, prbs1, prbs); % trl latencies
    rts(2,:) = cellfun(@(p1,p2) {nr(1, tlat>p1 & tlat<=p2)}, prbs1, prbs); % rs
    rts(3,:) = cellfun(@(p1,p2) {nr(2, tlat>p1 & tlat<=p2)}, prbs1, prbs); % rt


    % linear interpolation of missing responses
    rts(4,:) = cellfun(@(p1,p2) {interp1(p1(~isnan(p2)), p2(~isnan(p2)), p1, 'linear')}, rts(1,:), rts(2,:) );
    rts(5,:) = cellfun(@(p1,p2) {interp1(p1(~isnan(p2)), p2(~isnan(p2)), p1, 'linear')}, rts(1,:), rts(3,:) );


    % collate trial table
    tt = [];
    for px = 1:size(rts,2)  
        tt = [tt; [ repmat([snum, px, probe_res(px,[1,4,5,32,38])], numel(rts{1,px}), 1), ...
            rts{1,px}'/fs, (rts{1,px}-prbs{px})'/fs, reshape([rts{[2,4,3,5],px}], numel(rts{1,px}), []) ] ];
    end
    trl_tab = [trl_tab; [tt, [1:size(tt,1)]', test_res(:, [4,5,11,12]) ] ]; 


    % epoch-level summaries of mean and cv for interpolated RS and IBI series
    ep_10 = [prbs{:}]-fs*10;
    ep_30 = [prbs{:}]-fs*30;
    ep_off = [prbs{:}];
    mucv = nan(numel(prbs),12);
    for ex = 1:size(mucv,1)
        
        mucv(ex,1:2) = [snum, ex];
        
        % rs -- means on real responses, not interp
        mucv(ex,3) = mean( rts{2,ex}(rts{1,ex}>ep_30(ex) & rts{1,ex}<ep_off(ex)), 'omitnan' ); % ep_off should be redundant here but will be handy for ibi
        mucv(ex,4) = (std( rts{4,ex}(rts{1,ex}>ep_30(ex) & rts{1,ex}<ep_off(ex)), 'omitnan' ) / mucv(ex,3) );
        mucv(ex,5) = mean( rts{2,ex}(rts{1,ex}>ep_10(ex) & rts{1,ex}<ep_off(ex)), 'omitnan' ); % repeat for 10 s epochs
        mucv(ex,6) = (std( rts{4,ex}(rts{1,ex}>ep_10(ex) & rts{1,ex}<ep_off(ex)), 'omitnan' ) / mucv(ex,5) );

        % ibi
        mucv(ex,7) = mean( ibs{3,ex}(tmz{1,ex}(2):1/rs:tmz{1,ex}(end)>ep_30(ex)/fs & tmz{1,ex}(2):1/rs:tmz{1,ex}(end)<ep_off(ex)/fs) );
        mucv(ex,8) = (std( ibs{3,ex}(tmz{1,ex}(2):1/rs:tmz{1,ex}(end)>ep_30(ex)/fs & tmz{1,ex}(2):1/rs:tmz{1,ex}(end)<ep_off(ex)/fs) ) / mucv(ex,7) );
        mucv(ex,9) = mean( ibs{3,ex}(tmz{1,ex}(2):1/rs:tmz{1,ex}(end)>ep_10(ex)/fs & tmz{1,ex}(2):1/rs:tmz{1,ex}(end)<ep_off(ex)/fs) );
        mucv(ex,10) = (std( ibs{3,ex}(tmz{1,ex}(2):1/rs:tmz{1,ex}(end)>ep_10(ex)/fs & tmz{1,ex}(2):1/rs:tmz{1,ex}(end)<ep_off(ex)/fs) ) / mucv(ex,9) );
        
        mucv(ex,11) = mean( (ibs{3,ex}( tmz{1,ex}(2):1/rs:tmz{1,ex}(end)>ep_30(ex)/fs & tmz{1,ex}(2):1/rs:tmz{1,ex}(end)<ep_off(ex)/fs)-muIBI)/sdIBI );
        mucv(ex,12) = mean( (ibs{3,ex}( tmz{1,ex}(2):1/rs:tmz{1,ex}(end)>ep_10(ex)/fs & tmz{1,ex}(2):1/rs:tmz{1,ex}(end)<ep_off(ex)/fs)-muIBI)/sdIBI );

    end
    interpt = [interpt; mucv ];
   

    % output cells and aux variables
    save(['..' filesep 'ibi' filesep ['ibi_trl_' sname] ], 'ibs', 'rts', 'tmz', 'fs', 'rs', 'dfs', 'rt_thresh', 'prbs', 'tlat');

end

% label & convert IBI table
ibi_labs = {'subj_id', 'probe_num', 'probe_nested', 'block_num', 'stim_type', 'state', 'vigil', 'latency', 'time_to_probe', 'ibi', 'ibi_z'};
ibi_tab = array2table(ibi_tab, 'VariableNames', ibi_labs);

% label & convert trial/RT table
trl_labs = {'subj_id', 'probe_num', 'probe_nested', 'block_num', 'stim_type', 'state', 'vigil', 'latency', 'time_to_probe', ...
    'rspeed', 'rsInterp', 'rtime', 'rtInterp', 'trial_num', 'trial_nested', 'stim_num', 'trial_type', 'resp_type' };
trl_tab = array2table(trl_tab, 'VariableNames', trl_labs);

% deduce stim, trial and response (SDT) types
ibi_tab.stim_type(ibi_tab.stim_type==1)='F';
ibi_tab.stim_type(ibi_tab.stim_type==2)='D';
ibi_tab.stim_type = char(ibi_tab.stim_type);

trl_tab.stim_type(trl_tab.stim_type==1)='F';
trl_tab.stim_type(trl_tab.stim_type==2)='D';
trl_tab.stim_type = char(trl_tab.stim_type);

trl_tab.resp_type(trl_tab.trial_type==1)='C'; % nogo correct
trl_tab.resp_type(trl_tab.trial_type==0)='F'; % nogo error

trl_tab.resp_type(trl_tab.resp_type==0)='M'; % go error
trl_tab.resp_type(trl_tab.resp_type==1)='H'; % go correct
trl_tab.resp_type = char(trl_tab.resp_type);

trl_tab.trial_type(ismember(trl_tab.resp_type, ['H', 'M'])) = 'G';
trl_tab.trial_type(ismember(trl_tab.resp_type, ['F', 'C'])) = 'N';
trl_tab.trial_type = char(trl_tab.trial_type);

% output tables
writetable(ibi_tab, 'WIM_HB_ibi_data.txt')
writetable(trl_tab, 'WIM_HB_trl_data.txt')

writetable(array2table(interpt, ...
    'VariableNames', {'subj_id', 'probe_num', 'RSmu_30', 'RScv_30', 'RSmu_10', 'RScv_10',  'IBmu_30', 'IBcv_30',  'IBmu_10', 'IBcv_10', 'IBzu_30', 'IBzu_10' } ), ...
    'WIM_HB_IBI_RT_interp.txt')


