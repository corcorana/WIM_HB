function ECG = updateEvents(ECG, ann)

qrs = cell2mat(ann(1,:));

k = length(ECG.event);
evt = num2cell(qrs);
types = repmat({'QRS'},1, length(qrs));
urevents = num2cell(k+1:k+length(qrs));

[ECG.event(1,k+1:k+length(qrs)).latency] = evt{:};  % assign latencies
[ECG.event(1,k+1:k+length(qrs)).type] = types{:};   % assign types
[ECG.event(1,k+1:k+length(qrs)).urevent] = urevents{:}; % assign event index
[ECG.event(1,k+1:k+length(qrs)).duration] = deal(.5);

end