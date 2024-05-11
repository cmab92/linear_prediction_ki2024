function acf = estimate_acf(ts, L)
    %% in: (time series), (hypothesized model order)
    %% out: (acf estimate)
    %%
    D = L;
    N = D+L-1;
    X = hankel(ts(1:D), ts(D:end));
    v = 1./(max(max((abs((0:(N-1))-N/2)+D-N/2)/N, (D-L+1)/N), 0)');
    acf = v.*ifft(sum(fft(X, N).*conj(fft(X, N)), 2));
    acf = acf(1:L)/acf(1);
end

