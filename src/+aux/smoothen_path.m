%% path continuation - aux.smoothen_path
%  Interpolates data to archieve a smoothened path. 
%
%
%   Inputs:
%       X_data    -- contains a vector corresponding values or a 
%                    matrix which rows are vectors of corresponding values 
%       s_data    -- contains sample points for interpolation 
%
%   Optional Inputs: (call by value pairs)
%       'method'     -- interpolation method: e.g. linear, pchip, makima,
%                       spline
%       'number'     -- number of interpolation points. If passed,
%                       increment cannot be also an input! s_interp is
%                       calculated by linspace.
%       'increment'  -- increment between interpolation points If passed,
%                       number cannot be also an input! s_interp is
%                       calculated by s_start:increment:s_end.
%
%       --> By default 'method' is 'spline' and 'increment' is 1e-3.
%
%   Outputs:
%       X_interp     -- interpolated X_data
%       s_interp     -- interpolated s_data
%
%
%  See <a href="matlab:doc('interp1')">interp1</a> for detailed explanation of interpolation.
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>.
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   29.01.2022 - Tido Kubatschek
% 
function [X_interp,s_interp] = smoothen_path(X_data, s_data, varargin)
    %
    var_len = nargin - 2;
    %
    if mod(var_len,2) ~= 0
        error('There must be an even number of optional inputs!');
    end
    %
    %% default options
    %
    number = NaN;
    increment = 1e-3;
    method = 'spline';
    %
    %% check optional input
    %
    input_flag = 0;
    %
    for k = 1:2:var_len
        if strcmp(varargin{k}, 'number')
            if ~(isnumeric(varargin{k+1}) && varargin{k+1} == round(varargin{k+1}) && varargin{k+1} > 0)
                error('Wrong input for optional input number!');
            elseif input_flag ~= 0
                error('You cannot specify number and increment!');
            else
                number = varargin{k+1};
                increment = NaN;
                input_flag = 1;
            end
        elseif strcmp(varargin{k}, 'increment')
            if ~(isnumeric(varargin{k+1}) && varargin{k+1} > 0)
                error('Wrong input for optional input increment!');
            elseif input_flag ~= 0
                error('You cannot specify number and increment!');
            else
                increment = varargin{k+1};
                number = NaN;
                input_flag = 1;
            end
        elseif strcmp(varargin{k}, 'method')
            if ~ischar(varargin{k+1})
                error('Wrong input for optional input method!');
            else
                method = varargin{k+1};
                if strcmp(method, 'cubic')
                    warning('cubic doesnt work as the spacing is not equidistant, method is changed to spline');
                    method = 'spline';
                end
            end
        else
            error('Wrong optional input!');
        end
    end
    %
    if ~isnan(number)
        s_interp = linspace(s_data(1), s_data(end), number);
    else
        s_interp = s_data(1):increment:s_data(end);
    end
    %
    X_interp = (interp1(s_data, X_data.', s_interp, method)).';
end