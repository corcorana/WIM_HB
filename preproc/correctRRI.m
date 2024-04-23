function ECG = correctRRI(ECG, ann)

ECG.event = ECG.event(~strcmp({ECG.event.type}, 'ECG'));

qrs = cell2mat(ann(1,:));
k = length(ECG.event); % number of existing events in the EEG structure
urevents = num2cell(k+1:k+length(qrs));
evt = num2cell(qrs);
types = repmat({'ECG'},1, length(evt));

[ECG.event(1,k+1:k+length(qrs)).latency] = evt{:}; % assign latencies
[ECG.event(1,k+1:k+length(qrs)).type] = types{:}; % assign types
[ECG.event(1,k+1:k+length(qrs)).urevent] = urevents{:}; % assign event index
[ECG.event(1,k+1:k+length(qrs)).duration] = deal(1);

[~, idx] = sort([ECG.event.latency]);
ECG.event = ECG.event(idx);


end