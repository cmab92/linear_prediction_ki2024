function [pred, ts_complete] = predict_ts(ts, coeffs, n_steps)
    %% in: (known time series), (model coefficients), (# of time steps)
    %% out: (predicted values), (complete time series)
    %%
    ts = ts(:); coeffs = coeffs(:);
    order = length(coeffs);
    pred = zeros([n_steps, 1]);
    ts_complete = [ts; pred];
    for t = length(ts):length(ts)+n_steps-1
        ts_complete(t+1) = ts_complete(t-order+1:t)'*coeffs;
    end
    pred = ts_complete(end-n_steps+1:end);
end