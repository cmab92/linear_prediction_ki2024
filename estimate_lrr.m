function coeffs = estimate_lrr(eigvecs)
    %% in: (eigenvectors of covariance matrix)
    %% out: (estimated linear recurrence relations)
    %%
    v = power(norm(eigvecs(end, :)), 2);
    coeffs = 1/(1-v)*(eigvecs(1:end-1, :))*eigvecs(end, :)';
    if sum(abs(roots(coeffs))<1)>0
        disp("LRR instable.")
    end
end