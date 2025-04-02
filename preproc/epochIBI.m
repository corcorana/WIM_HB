function ibi = epochIBI(ECG)

% divide ECG into segements between successive probes
tcodes =  {'S 65' 'S 66' 'S 67' 'S 68' 'S 69' 'S 70' 'S 71' 'S 72' 'S 73', 'T  1'}; % triggers for trial onsets

% get latencies for events
trls = [ECG.event(ismember({ECG.event.type}, tcodes)).latency]; % trial onsets
prbs = [ECG.event(ismember({ECG.event.type}, {'S  3' 'P  1'})).latency]; % probe onsets
rpks = [ECG.event(ismember({ECG.event.type}, 'QRS')).latency]; % R-peak events

% sort trials according to probe intervals
trls = trls( ~(trls>prbs(end)) ); % exclude trial latencies after final probe
pmat = trls>prbs'; % identify trials following nth probe
csum = trls(logical(sum(cumsum(pmat,2)==1))); % latency of first trial following nth probe

% sample R-peak event preceding 1st trial to follow previous probe & 1st R-peak following next probe 
% (i.e., 1st IBI contains 1st trial; last IBI contains probe)
ibi = cell(1, length(prbs));
fprb = rpks < prbs(1); 
ibi{1,1} = rpks(rpks>=ECG.event(find(ismember({ECG.event.type}, tcodes),1, 'first')-1).latency & [true, fprb(1:end-1)]);
for ex = 1:length(csum) 
    fprb = rpks < prbs(ex+1);
    ibi{1,ex+1} = rpks(rpks>=ECG.event( find([ECG.event.latency]==csum(ex) & ...
        ismember({ECG.event.type}, tcodes))-1).latency & [true, fprb(1:end-1)]);
end

end