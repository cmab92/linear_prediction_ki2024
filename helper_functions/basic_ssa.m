function x_rec = basic_ssa(x, L, idx)
    %% in: (time series), (embedding dimension), (index of eigentriples)
    %% out: (reconstructed singal)
    % (almost) vanilla SSA
    % embedding:
    X = hankel(x(1:L), x(L:end));
    [m, n] = size(X);
    % decomposition / SVD
    [U, S, V] = svd(X);
    X = U(:, idx)*S(idx, idx)*V(:, idx)';
    % diagonal averaging (or rather antidiagonal averaging)
    % Y = diagonal_averaging(X);
    % x_rec = recover(Y)';
    X = fliplr(X);
    diag_index = -m+1:n-1;
    x_rec = zeros([m+n-1, 1]);
    for i = 1:m+n-1
        x_rec(i) = mean(diag(X, diag_index(i)));
    end
    % vanilla: x_rec = flip(x_rec), X_avg = hankel(x_rec(1:L), x_rec(L:end))
    % recontruction:
    x_rec = flip(x_rec);
    % vanilla: x_rec = [X_avg(1, :) X_avg(2:end, end)'];
end