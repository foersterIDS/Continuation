%% path continuation - aux.smoothenPath
%  Interpolates data to archieve a smoothened path. 
%
%
%   Inputs:
%       XData    -- contains a vector corresponding values or a 
%                    matrix which rows are vectors of corresponding values 
%       sData    -- contains sample points for interpolation 
%
%   Optional Inputs: (call by value pairs)
%       'method'     -- interpolation method: e.g. linear, pchip, makima,
%                       spline
%       'number'     -- number of interpolation points. If passed,
%                       increment cannot be also an input! sInterp is
%                       calculated by linspace.
%       'increment'  -- increment between interpolation points If passed,
%                       number cannot be also an input! sInterp is
%                       calculated by sStart:increment:sEnd.
%
%       --> By default 'method' is 'spline' and 'increment' is 1e-3.
%
%   Outputs:
%       XInterp     -- interpolated XData
%       sInterp     -- interpolated sData
%
%
%  See <a href="matlab:doc('interp1')">interp1</a> for detailed explanation of interpolation.
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>.
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   29.01.2022 - Tido Kubatschek
% 
function [XInterp,sInterp] = smoothenPath(XData, sData, varargin)
    %
    varLen = nargin - 2;
    %
    if mod(varLen,2) ~= 0
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
    inputFlag = 0;
    %
    for k = 1:2:varLen
        if strcmp(varargin{k}, 'number')
            if ~(isnumeric(varargin{k+1}) && varargin{k+1} == round(varargin{k+1}) && varargin{k+1} > 0)
                error('Wrong input for optional input number!');
            elseif inputFlag ~= 0
                error('You cannot specify number and increment!');
            else
                number = varargin{k+1};
                increment = NaN;
                inputFlag = 1;
            end
        elseif strcmp(varargin{k}, 'increment')
            if ~(isnumeric(varargin{k+1}) && varargin{k+1} > 0)
                error('Wrong input for optional input increment!');
            elseif inputFlag ~= 0
                error('You cannot specify number and increment!');
            else
                increment = varargin{k+1};
                number = NaN;
                inputFlag = 1;
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
    
    %% Iteration der Interpolation
    diffS = inf;
    tolS = 1e-6;
    maxIt = 10;
    kIt = 0;

    while diffS > tolS && kIt < maxIt
        if ~isnan(number)
            sInterp = linspace(sData(1), sData(end), number);
        else
            sInterp = sData(1):increment:sData(end);
        end

        % Calculate interpolated points of curve
        XInterp = (interp1(sData, XData.', sInterp, method)).';

        % With those recalculate arclength
        sInterp = [0,cumsum(sqrt(sum((XInterp(:,2:end) - XInterp(:,1:end-1)).^2,1)))];

        % caculate difference in covered arclength
        diffS = abs(sData(end) - sInterp(end));

        % update
        sData = sInterp;
        XData = XInterp;

        kIt = kIt + 1;
    end
end