%% path continuation - numeric_jacobian
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%   16.09.2020 - Tido Kubatschek 
%
function jac = numeric_jacobian(f, x, varargin)
    % Calculate Jacobian of function f at given x
    %
    % varargin consists of
    %   'fx' followed by optional f(x)
    %   'is' followed by starting value for computation of the jacobian
    
    epsilon = 1e-6;
    epsilon_inv = 1/epsilon;
    nx = length(x); % Dimension of the input x;
    
    % test for variabel inputs
    vl = numel(varargin);
    if vl > 6
        error('Too many inputs!')
    end
    
    fi = 0;
    ii = false;
    diff_quo = 0;
    
    for ki = 1:vl
        if strcmp(varargin{ki},'central_value')
            f0 = varargin{ki+1};
            fi = 1;
        end
        if strcmp(varargin{ki}, 'derivative_dimensions')
            is = varargin{ki+1};
            ii = true;
        end
        if strcmp(varargin{ki}, 'diffquot')
            diff_str = varargin{ki+1};
            if diff_str.central
                diff_quo = 1;    
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
        x_ = x;
        x_(is(i)) =  x(is(i)) + epsilon;
        x_r = x;
        x_r(is(i)) = x(is(i)) - epsilon;
        if diff_quo == 0
            %% forward difference quotient
            jac(:, i) = (feval(f, x_) - f0) .* epsilon_inv;
        else
            %% central difference quotient
            jac(:, i) = (feval(f,x_) - feval(f,x_r)) .* (0.5* epsilon_inv);
        end
    end
    
end