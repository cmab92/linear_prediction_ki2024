function z = linconv_projection(ts, eigvecs)
    %% in: (time series), (eigenvectors of (auto)covariance matrix)
    %% out: (projected time series)
    %%
    L = size(eigvecs, 1);
    n_vecs = size(eigvecs, 2);
    D = length(ts);
    N = D+L-1;
    %
    z = zeros(size(ts));
    ts_fd = fft(ts, N);
    eigvecs_fd = fft(eigvecs, N);
    for i = 1:n_vecs
        temp = ifft(eigvecs_fd(:, i).*ts_fd);                  % circ. conv
        temp = real(ifft(conj(eigvecs_fd(:, i)).*fft(temp)));   % backward
        z = z + temp(1:D);
    end
    z = z/L;
end