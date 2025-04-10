function t_perms = lmeEEG_permute(mEEG, X, tab, nperms)

rperms = lmeEEG_permutations(tab.sid, nperms); % nperms within-subjects permutations of X (for stimuli-within-condition designs) 
t_perms = nan(nperms, size(mEEG,1), size(mEEG,2), size(X,2)); % Initialize t-map

for np = 1:nperms
    fprintf('Permutation %g/%g ...\n', np, nperms)
    XX = X(rperms(:,np),:,:);
    for ch = 1:size(mEEG,1) % channel
        % exclude missing entries in design matrix (rejected chans / absent HEPs)
        pexc = isnan(XX(:,1,ch));
        Xx = XX(~pexc,:,ch);
        % ... and corresponding mEEG
        excl = isnan(X(:,1,ch));
        mEEGx = mEEG(:,:,~excl);
        parfor sp = 1:size(mEEG,2) % sample point (time or freq)
            EEGx = squeeze(mEEGx(ch,sp,:));
            t_perms(np,ch,sp,:) = lmeEEG_regress(EEGx, Xx);
        end
    end
end


end