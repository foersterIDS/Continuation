%% path continuation - aux.numericJacobian
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%   16.09.2020 - Tido Kubatschek 
%
function [jac,f0] = numericJacobian(f, x, varargin)
    % Calculate Jacobian of function f at given x
    %
    % varargin consists of
    %   'fx' followed by optional f(x)
    %   'is' followed by starting value for computation of the jacobian
    
    epsilon = 1e-6;
    epsilonInv = 1/epsilon;
    nx = length(x); % Dimension of the input x;
    
    % test for variabel inputs
    vl = numel(varargin);
    if vl > 6
        error('Too many inputs!')
    end
    
    fi = 0;
    ii = false;
    diffQuo = 0;
    
    for ki = 1:vl
        if strcmp(varargin{ki},'centralValue')
            f0 = varargin{ki+1};
            fi = 1;
        end
        if strcmp(varargin{ki}, 'derivativeDimensions')
            is = varargin{ki+1};
            ii = true;
        end
        if strcmp(varargin{ki}, 'diffquot')
            diffStr = varargin{ki+1};
            if diffStr.central
                diffQuo = 1;    
            end
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
            jac(:, i) = (feval(f, xn) - f0) .* epsilonInv;
        else
            %% central difference quotient
            jac(:, i) = (feval(f,xn) - feval(f,xr)) .* (0.5* epsilonInv);
        end
    end
    
end