%% path continuation - aux.numericJacobian
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%   16.09.2020 - Tido Kubatschek 
%
function [jac,f0] = numericJacobian(f, x, NameValueArgs)
    % Calculate Jacobian of function f at given x
    %
    % varargin consists of
    %   'fx' followed by optional f(x)
    %   'is' followed by starting value for computation of the jacobian

    arguments
        f (1,1) function_handle
        x (:,1) double
        NameValueArgs.centralValue (:,1) double
        NameValueArgs.derivativeDimensions (1,:) double {mustBeInteger,mustBeGreaterThanOrEqual(NameValueArgs.derivativeDimensions,1)}
        NameValueArgs.diffQuot (1,1) struct
    end
    
    epsilon = 1e-6;
    epsilonInv = 1/epsilon;
    nx = length(x); % Dimension of the input x;
    
    fi = 0;
    ii = false;
    diffQuo = 0;

    if isfield(NameValueArgs,'centralValue')
        f0 = NameValueArgs.centralValue;
        fi = 1;
    end
    if isfield(NameValueArgs,'derivativeDimensions')
        is = NameValueArgs.derivativeDimensions;
        ii = true;
    end
    if isfield(NameValueArgs,'diffQuot')
        diffStr = NameValueArgs.diffQuot;
        if diffStr.central
            diffQuo = 1;    
        end
    end
    
    if fi == 0
        f0 = feval(f, x); % caclulate f0, when no perturbation happens
    end
    
    if ~ii
        is = 1:nx;
    end
        
    % Do perturbation
    jac = NaN(length(f0),length(is));
    for i = 1:length(is)
        xn = x;
        xn(is(i)) =  x(is(i)) + epsilon;
        xr = x;
        xr(is(i)) = x(is(i)) - epsilon;
        if diffQuo == 0
            %% forward difference quotient
            jac(:, i) = (f(xn) - f0) .* epsilonInv;
        else
            %% central difference quotient
            jac(:, i) = (f(xn) - f(xr)) .* (0.5* epsilonInv);
        end
    end
    
end