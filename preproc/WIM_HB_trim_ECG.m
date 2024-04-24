%%
clear
close all

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_eeglab))

%% batch process
subjs = dlmread(['..' filesep 'subjs.txt']);


for ix = 1 :length(subjs)

    snum = subjs(ix);
    sname = ['MWI' num2str(snum)];
    fprintf(['Loading ' sname '...\n'])

    % load corrected ECG data
    ECG = pop_loadset( 'filename', [sname, '_ibi_corrected.set'], 'filepath', ['..' filesep 'ibi']);
    
    % downsample   
    ECG = pop_resample(ECG, 250);

    % trim data per EEG preprocessing
    bsta = find(ismember({ECG.event.type}, 'B  1'));
    bsta = bsta(2:end);
    bend = find(ismember({ECG.event.type}, 'K  1'));
    bend = bend(1:end-1);
    for bx = length(bend):-1:1
        if (ECG.event(bsta(bx)).latency) - (ECG.event(bend(bx)).latency+ECG.srate*2) > ECG.srate
            ECG = pop_select( ECG, 'nopoint', [ECG.event(bend(bx)).latency+ECG.srate*2, ECG.event(bsta(bx)).latency-1] );
        end
    end
    psta = find(contains({ECG.event.type}, 'P  '));
    pend = find(matches({ECG.event.type}, 'C  1'));
    for px = length(psta):-1:1
    	ECG = pop_select( ECG, 'nopoint',[ECG.event(psta(px)).latency+ECG.srate*2, ECG.event(pend(px)).latency-1] );
    end

    % save trimmed file
    ECG = pop_saveset( ECG, 'filename', [sname, '_ibi_corrected_trimmed.set'], 'filepath', ['..' filesep 'ibi']);

end
