function x_new = MaxVariance_function(obj, A, strategy)
%MAXVARIANCE_FUNCTION Select the point with maximum Kriging prediction variance.
%
% This adaptive sampling function selects a new sample point by maximizing
% the Kriging prediction variance.
%
% Inputs:
%   obj      - Current OK_model object
%   A        - Sampling domain
%              A(1,:) = lower bounds
%              A(2,:) = upper bounds
%   strategy - Optimization strategy, for example 'MS'
%
% Output:
%   x_new    - New sample point in unity space

    % Lower and upper bounds
    lb = A(1,:);
    ub = A(2,:);

    % Linear constraints: none
    AA = [];
    b = [];
    Aeq = [];
    beq = [];

    % Nonlinear constraints: none
    nonlcon = [];

    % Objective function:
    % optimizationTools minimizes the objective,
    % so we minimize negative prediction variance.
    fun = @(x) MaxVariance_objective(obj, x);

    % Run optimization
    x_new = optimizationTools(fun, strategy, AA, b, Aeq, beq, lb, ub, nonlcon);

end


function f = MaxVariance_objective(obj, x)
%MAXVARIANCE_OBJECTIVE Objective function for maximum variance sampling.
%
% The second output of obj.predict(x) is assumed to be the Kriging
% prediction variance or MSE.
%
% Since the optimizer minimizes f, we use:
%
%   f = -variance
%
% which is equivalent to maximizing the prediction variance.

    % Kriging prediction
    [~, variance] = obj.predict(x);

    % Avoid numerical negative variance
    variance = max(variance, 0);

    % Minimize negative variance = maximize variance
    f = -variance;

end