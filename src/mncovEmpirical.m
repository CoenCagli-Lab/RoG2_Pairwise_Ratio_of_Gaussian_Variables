function [muEmp, sigmaEmp, corrEmp] = mncovEmpirical(spikes)
%% mncovEmpirical Computes the empirical mean and covariance matrices for pairwise neural response
%   INPUT
%       spikes (number of pairs, 2, number of stimuli, number of trials) - pairwise neural responses
%   OUTPUT
%       muEmp (number of pairs, 2, number of stimuli) - empirical mean
%       sigmaEmp (number of pairs, number of stimuli, 2, 2) - empirical covariance matrix
%       corrEmp (number of pairs, number of stimuli) - empirical noise correlation
%   NOTE: If there is only 1 pair, the first dimension will be removed
%%
    if numel(size(spikes)) == 4
        [nPairs, ~, nStims, ~] = size(spikes);
        muEmp = NaN(nPairs, 2, nStims);
        sigmaEmp = NaN(nPairs, nStims, 2, 2);
        corrEmp = NaN(nPairs,nStims);
        for iPair = 1:nPairs
            [muEmp(iPair, :, :), sigmaEmp(iPair, :, :, :)] = ...
                internal_mncov(squeeze(spikes(iPair, :, :, :)));
            corrEmp(iPair,:) = squeeze(sigmaEmp(iPair,:,1,2))./sqrt(squeeze(sigmaEmp(iPair,:,1,1)).*squeeze(sigmaEmp(iPair,:,2,2)));
        end

    elseif numel(size(spikes)) == 3
        [muEmp, sigmaEmp] = internal_mncov(spikes);
        corrEmp = squeeze(sigmaEmp(:,1,2))./sqrt(squeeze(sigmaEmp(:,1,1)).*squeeze(sigmaEmp(:,2,2)));
    else
        error('Wrong spike format');
    end


end

function [mu, sig] = internal_mncov(spikes)
    nStims = size(spikes, 2);
    mu = squeeze(mean(spikes, 3, 'omitnan'));
    sig = NaN(nStims, 2, 2);
    for iStim = 1:nStims
        sig(iStim, :, :) = cov(squeeze(spikes(:, iStim, :))', 'omitrows');
    end
end
