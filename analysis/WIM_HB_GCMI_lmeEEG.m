%%
clear
close

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(genpath(path_lmeEEG))
addpath(genpath(path_TFCE))

%% batch process

load('WIM_HB_GCMI_NMD.mat')

% concatenate data across subjects
DAT = cat(5, MI{:,1});
P = cat(1, MI{:,2});

% factor info
sid = categorical( P(:,1) );
pid = categorical( P(:,2) );
site = categorical(P(:,1)<100, 0:1, {'MBI', 'PBI'} );
stim = categorical(P(:,3), 1:2, {'FACE', 'DIGIT'} );
ms = categorical(P(:,4), 1:3, {'ON', 'MW', 'MB'}); 
vig = categorical(P(:,5), 1:4, {'V4', 'V3', 'V2', 'V1'}); % reverse code so V4 = 'Extremely Alert'

TAB = table(sid, pid, site, stim, ms, vig);

% analysis
mform = 'y~site+stim+x+(1|sid)+(1|pid)'; %  model formula 
nperms = 1000;

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
    dat = DAT(:,:,:,:,eps);
    
    for dx = 1:size(dat,1) % directionality (1=hb, 2=bh)
    
        d = squeeze(mean(dat(dx,:,:,:,:),4));
    
        % fit mixed models & perform mass univariate regression
        [Results(sx,dx).t_obs, Results(sx,dx).betas, Results(sx,dx).se, Results(sx).cnams, mEEG, X] = lmeEEG_fitMods(d, tab, mform);
        
        % permute design matrix
        Results(sx,dx).t_perms = lmeEEG_permute(mEEG, X, tab, nperms);

        % apply TFCE
        for jx = 2:size(X,2)
            Results(sx,dx).(matlab.lang.makeValidName(Results(sx).cnams{jx})) = lmeEEG_TFCE( ...
                squeeze(Results(sx,dx).t_obs(:,:,jx)), squeeze(Results(sx,dx).t_perms(:,:,:,jx)), chanlocs, [0.66 2]);
        end

    end

end

% save results
save('WIM_HB_GCMI_lmeEEG.mat', 'Results', 'TAB', 'mform', 'fbands', 'cfrex', 'chanlocs')
