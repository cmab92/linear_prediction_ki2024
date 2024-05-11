function coeffs = estimate_wiener(x, y, order)
    %% in: (input vector), (output vector), (length of filter)
    %% out: (adaptive filter coefficients)
    %%
    X = hankel(x(1:order), x(order:end));
    coeffs = y'*pinv(X);
    if sum(abs(roots(coeffs))<1)>0
        disp("Wiener filter instable.")
    end
end