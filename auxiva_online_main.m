clear;

resampFreq = 16e3;

[srcSig(:,:,1), sampFreq] = audioread('drums.wav'); % signal x channel x source (source image)
[srcSig(:,:,2), sampFreq] = audioread('piano.wav'); % signal x channel x source (source image)
srcSigResample(:,:,1) = resample(srcSig(:,:,1), resampFreq, sampFreq, 100); % resampling for reducing computational cost
srcSigResample(:,:,2) = resample(srcSig(:,:,2), resampFreq, sampFreq, 100); % resampling for reducing computational cost

% Mix source images of each channel to produce observed mixture signal
mixSig(1, :) = srcSigResample(:,1,1) + srcSigResample(:,1,2);
mixSig(2, :) = srcSigResample(:,2,1) + srcSigResample(:,2,2);
if abs(max(max(mixSig))) > 1 % check clipping
    error('Cliping detected while mixing.\n');
end

[nCh, nSamples] = size(mixSig);

nFft = 512;
nOverlap = 2;
nDelay = nOverlap - 1;
nHop = nFft / nOverlap;
% winFunc = hamming(nFft, "periodic").';
w=hann(nFft, 'periodic');
winFunc=sqrt(w*2.0*nHop/nFft).';

nFrame = floor(nSamples / nHop) - 1;

sigDelay = zeros(nHop*nDelay, nCh).';
sigOut = zeros(nSamples, nCh).';

processedSig = zeros(nFft, nCh).';
processedSigDelay = zeros(nHop*nDelay, nCh).';

sigPartial = zeros(nFft, nCh);

scalingtmp = zeros(nFft*2, 1).';
scalingtmp(1:nFft) = winFunc.^2;
scalingtmp(nHop+1:nHop+nFft) = scalingtmp(nHop+1:nHop+nFft) + winFunc.^2;

scaling = 1./scalingtmp(nHop+1:nHop+nHop);
nBins = nFft/2+1;


% Initialization
[n_frame, n_freq, n_src] = size(mix);
eye_matrix = eye(n_src);

W = repmat(eye_matrix, [n_frame, n_freq, 1, 1]);
cov = repmat(eye_matrix, [n_frame, n_src, n_freq, 1, 1]) * cov_scale;
est = zeros(size(mix));
est_out = zeros(size(mix));

cont = containers.Map({'Gauss', 'Laplace'}, ...
    { @(y) sum(abs(y).^2, 2) / n_freq, @(y) 2 * sum(abs(y), 2) });

% Iteration
for t = n_blocks:n_frame
    W(t, :, :, :) = W(t - 1, :, :, :);

    t_block = t - n_blocks + 1:t;

    for k = 1:n_iter
        % Update source model
        est(t_block, :, :, :) = W(t_block, :, :, :) * permute(mix(t_block, :, :), [1, 3, 2]);

        r = cont(update_source_model(est(t_block)));
        for s = 1:n_src
            % Update source model
            est(t_block, :, :, :) = W(t_block, :, s, :) * permute(mix(t_block, :, :), [1, 3, 2]);
            r = cont(update_source_model(est(t_block)));

            % Update weighted covariance
            cov_new = mean(bsxfun(@times, mix(t_block, :, :), 1 ./ r) * permute(conj(mix(t_block, :, :)), [1, 3, 2]), 1);
            cov(t, s, :, :, :) = forget_param * cov(t - n_blocks, s, :, :, :) + (1 - forget_param) * cov_new;

            % Update demixing vector
            W(t, :, :, :) = update_spatial_model(cov(t, :, :, :, :), W(t, :, :, :), 'row_idx', s, 'method', update_demix_filter);
        end
    end

    W_pb = projection_back_frame(W(t, :, :), ref_mic);
    est_out(t, :, :) = demix(mix(t, :, :), W_pb);
end

function W_proj = projection_back_frame(W, ref_mic)
    [n_freq, ~, ~] = size(W);
    W_proj = W;

    A = inv(W);
    for f = 1:n_freq
        eA = diag(A(f, ref_mic, :));
        W_proj(f, :, :) = eA * W(f, :, :);
    end
end

