% Auxiliary-function-based independent vector analysis
function [W, est_out] = auxiva_online(mix, update_demix_filter, update_source_model, n_blocks, forget_param, n_iter, ref_mic, cov_scale)
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

