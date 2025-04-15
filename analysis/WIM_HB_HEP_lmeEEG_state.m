%%
clear
close

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_lmeEEG))
addpath(genpath(path_TFCE))

%% batch process

load('WIM_HB_HEP.mat')

% concatenate data across subjects
DAT = cat(3, allHEP{:,1});
P = cat(1, allHEP{:,2});

% factors
sid = categorical( P(:,1) ); 
pid = P(:,2);
site = categorical(P(:,1)<100, 0:1, {'MBI', 'PBI'} );
stim = categorical(P(:,3), 1:2, {'FACE', 'DIGIT'} );
ms = categorical(P(:,4), 1:3, {'ON', 'MW', 'MB'});
vig = categorical(P(:,5), 1:4, {'V4', 'V3', 'V2', 'V1'}); % reverse code so V4 = 'Extremely Alert'

ibs = readtable("WIM_HB_IBI_pup_behav.csv");
TAB = [ table(sid, pid, site, stim, ms, vig), table(ibs.zuIBI, normalize(log(ibs.cvIBI)), 'VariableNames', {'muIBI', 'cvIBI'}) ];

% analysis 
mform = 'y~site+stim+muIBI+cvIBI+x+(1|sid)'; %  model formula 
nperms = 1000; % number of design matrix permutations
twin = t>=.2 & t<=.6; % time-window for TFCE (s)

Results = struct([]);

for sx = 1:2 % state var (1 = MS; 2 = VIG)

    % define state variable & exclude epochs with missing data
    if sx == 1
        u = unique(TAB.sid(TAB.ms~="ON")); % identify subjects that report OFF states
        eps = ismember(TAB.sid, u) & ~isundefined(TAB.ms); % exclude subjects that don't report OFF; any missing probe reports
        x = TAB.ms(eps);
    elseif sx == 2
        eps = ~isundefined(TAB.vig);
        x = TAB.vig(eps); % exclude missing VIG ratings
    end
    tab = [TAB(eps,:), table(x)];
    dat = DAT(:,:,eps);

    % fit mixed models & perform mass univariate regression
    [Results(sx).t_obs, Results(sx).betas, Results(sx).se, mEEG, X, mod] = lmeEEG_fitMods(dat, tab, mform);
    
    % permute design matrix
    Results(sx).t_perms = lmeEEG_permute(mEEG, X, tab, nperms);
 
    % apply TFCE
    for jx = 2:size(X,2)
        Results(sx).(matlab.lang.makeValidName(mod.CoefficientNames{jx})) = lmeEEG_TFCE( ...
            squeeze(Results(sx).t_obs(:,twin,jx)), squeeze(Results(sx).t_perms(:,:,twin,jx)), chanlocs, [0.66 2]);
    end

end

% save results
save('WIM_HB_HEP_lmeEEG.mat', 'Results', 'TAB', 'mform', 'twin', 't', 'chanlocs', '-v7.3')
