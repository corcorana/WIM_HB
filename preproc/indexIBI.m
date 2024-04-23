function ibi = indexIBI(ECG)

% divide ECG into segements between successive probes
t = regexp({ECG.event.type},'T  (\d)([FD])', 'tokens'); % don't mistake T  1 markers as trial markers
trls = [ECG.event(~cellfun(@isempty, t)).latency]; % trial events
prbs = [ECG.event(contains({ECG.event.type}, 'P')).latency]; % probe events
rpks = [ECG.event(ismember({ECG.event.type}, 'ECG')).latency]; % R-peak events
pmat = trls>prbs'; % trials after nth probe
csum = trls(logical(sum(cumsum(pmat,2)==1))); % latency of first trial following nth probe

% sample R-peak event preceding first trial to follow previous probe &
% first R-peak following next probe (1st IBI contains 1st trial; last IBI contains probe)
ibi = cell(1, length(prbs));
fprb = rpks < prbs(1); 
ibi{1,1} = rpks(rpks>=ECG.event(find(~cellfun(@isempty, t),1, 'first')-1).latency & [true, fprb(1:end-1)]);
for ex = 1:length(csum)
    fprb = rpks < prbs(ex+1);
    ibi{1,ex+1} = rpks(rpks>=ECG.event( find([ECG.event.latency]==csum(ex) & ...
        contains({ECG.event.type},'T  ')) -1).latency & [true, fprb(1:end-1)]);
end

end