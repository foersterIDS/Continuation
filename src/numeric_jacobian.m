%% path continuation - numeric_jacobian
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%
function jac = numeric_jacobian(f, x, fx)
    % Calculate Jacobian of function f at given x
    % fx is optional f(x)
    
    epsilon = 1e-6; 
    epsilon_inv = 1/epsilon;
    nx = length(x); % Dimension of the input x;
    if ~(nargin>2)
        f0 = feval(f, x); % caclulate f0, when no perturbation happens
    else
        f0 = fx;
    end
    
    % Do perturbation
    for i = 1 : nx
        x_ = x;
        x_(i) =  x(i) + epsilon;
        jac(:, i) = (feval(f, x_) - f0) .* epsilon_inv;
    end
    
end