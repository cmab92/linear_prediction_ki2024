function z = ssa_projection(ts, eigvecs)
    %% in: (time series), (eigenvectors of (auto)covariance matrix)
    %% out: (projected time series)
    %%
    L = size(eigvecs, 1);
    n_vecs = size(eigvecs, 2);
    D = length(ts);
    K = min([L, D-L+1]);
    w = K*ones([D, 1]); w(1:K) = (1:K); w(end-K+1:end) = (K:-1:1);
    w = 1./w;
    %
    z = zeros(size(ts));
    ts_fd = fft(ts, D);
    eigvecs_fd = (fft(flip(eigvecs), D));
    for i = 1:n_vecs
        temp = real(ifft(eigvecs_fd(:, i).*ts_fd));               % circ. conv
        temp(1:L-1) = 0;                                          % truncation (valid convolution) & padding for subsequent backward filtering
        temp = real(ifft(conj(eigvecs_fd(:, i)).*fft(temp)));     % backward
        z = z + temp;
    end
    z = w.*z;     % diag(w) is the inverse frame operator
end