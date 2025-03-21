function err = LogCosh(p, F_params, prices, T1, T2, m, n, X_t, Y_t)
% This function calculates the Log-Cosh Loss between the
% Future price values predicted by the Arithmetic model and observed
% market prices.
%
% Parameters:
% p:        Vector of model parameters
% F_params: Function handle for computing forward prices based on the 
%           Arithmetic model
% prices:   Observed market prices
% T1:       Vector of start delivery times
% T2:       Vector of end delivery times (settlements)
% m:        Number of Gaussian processes in the model.
% n:        Number of Jump processes in the model.
% X_t:      Current values of Gaussian processes.
% Y_t:      Current values of Jump processes.
%
% Returns:
% err       The Log-Cosh Loss between model-predicted and observed prices

% Calculate the Log-Cosh Loss by iterating over all Future contracts
err = mean(arrayfun(@(i) ...
    log(cosh(F_params(0, T1(i), T2(i), p, m, n, X_t, Y_t) ...
    - prices(i))), ...
    1:length(T1)));

end

