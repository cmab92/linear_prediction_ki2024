function coeffs = estimate_ar(acf, order)
    %% in: (empirical autocorrelation sequence), (hypothesized model order)
    %% out: (estimated AR model coefficients)
    %%
    coeffs = levinson(acf, order);
    coeffs = -flip(coeffs(2:end));  % cf. https://de.mathworks.com/help/signal/ref/levinson.html
    if sum(abs(roots(coeffs))<1)>0
        disp("AR polynomial instable.")
    end
end