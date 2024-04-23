%%
clear
close all

%% paths
run(['..' filesep 'localdef_WIM_HB'])

%% batch process
snum = dlmread(['..' filesep 'subjs.txt']);

sex = []; age = [];
for ix = 1 :length(snum)
    
    % load behav data (for probes)
    behav_name = dir([path_behav filesep '*s' num2str( snum(ix) ) '*.mat']);
    load([path_behav filesep behav_name.name]);

    sex = [sex; SubjectInfo.Gender];
    age = [age; str2double(SubjectInfo.Age)];

end

tab = table(snum, sex, age);

table( sum(tab.sex == 'F'), sum(tab.sex == 'M'), min(tab.age), max(tab.age), mean(tab.age), std(tab.age), ...
    'VariableNames', {'n female', 'n male', 'min age', 'max age', 'mean age', 'sd age'}) 


