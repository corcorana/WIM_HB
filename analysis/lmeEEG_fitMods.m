function [t_obs, betas, se, cnams, mEEG, X] = lmeEEG_fitMods(dat, tab, mform)

% Extract number of cols required for design matrix X
y = squeeze(dat(1,1,:));
df = [table(y), tab];
mod = fitlme(df, mform);
cnams = mod.CoefficientNames;

% Conduct mixed models on each channel/timepoint combination.     
mEEG = nan(size(dat)); 
X = nan(size(dat,3), length(cnams), size(dat,1));
for ch = 1:size(dat,1) % channel
    fprintf('Fitting mixed models [%g/%g]...\n', ch, size(dat,1))
    parfor sp = 1:size(dat,2) % sample point (time or freq)
        y = squeeze(dat(ch,sp,:));
        df = [table(y), tab];
        mod = fitlme(df, mform); 
        mEEG(ch,sp,:) = fitted(mod,'Conditional',0)+residuals(mod); % Extract marginal EEG
    end
    % Extract design matrix X for channel data (differs due to rejected chans)
    y = squeeze(dat(ch,1,:));
    df = [table(y), tab];
    mod = fitlme(df, mform);
    X(:,:,ch) = designMatrix(mod);
end

% Perform mass univariate linear regressions on “marginal” EEG data.
[t_obs, betas, se] = deal( nan(size(mEEG,1), size(mEEG,2), size(X,2)) );
for ch = 1:size(mEEG,1)
    % exclude missing entries in design matrix (rejected chans / epochs)
    excl = isnan(X(:,1,ch));
    Xx = X(~excl,:,ch);
    % ... and corresponding mEEG
    mEEGx = mEEG(:,:,~excl);
    fprintf('Obtaining statistics [%g/%g]...\n', ch, size(dat,1))
    parfor sp = 1:size(mEEGx,2)
        EEGx = squeeze(mEEGx(ch,sp,:));
        [t_obs(ch,sp,:), betas(ch,sp,:), se(ch,sp,:)] = lmeEEG_regress(EEGx, Xx);
    end
end

end