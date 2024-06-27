%% path continuation - aux.smoothenPath
%  Interpolates data to archieve a smoothened path.
%
%  See <a href="matlab:doc('interp1')">interp1</a> for detailed explanation of interpolation.
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>.
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   29.01.2022 - Tido Kubatschek
% 
function [varInterp,paraInterp,sInterp] = smoothenPath(varData, paraData, sData, NameValueArgs)
    arguments
        varData (:,:) double {mustBeReal}
        paraData (:,:) double {mustBeReal}
        sData (1,:) double {mustBeIncreasing(sData)} = []
        NameValueArgs.number double {mustBeInteger,mustBeGreaterThan(NameValueArgs.number,0),mustBeScalarOrEmpty} = []
        NameValueArgs.increment double {mustBeInteger,mustBeGreaterThan(NameValueArgs.increment,0),mustBeScalarOrEmpty} = []
        NameValueArgs.interpolationMethod (1,:) char {mustBeMember(NameValueArgs.interpolationMethod, {'linear','pchip','makima','spline'})} = 'spline'
        NameValueArgs.sTolerance (1,1) double {mustBeGreaterThan(NameValueArgs.sTolerance,0)} = 1e-2
        NameValueArgs.maxIterations (1,1) double {mustBeInteger,mustBeGreaterThan(NameValueArgs.maxIterations,0)} = 3
    end
    %% get optional inputs
    if ~isempty(NameValueArgs.number) && ~isempty(NameValueArgs.increment)
        warning('You cannot specify both NameValueArgs.number and NameValueArgs.increment. Input for NameValueArgs.increment is used.')
        numberOfInterpolationPoints = [];
        interpolationIncrement = NameValueArgs.increment;
    elseif ~isempty(NameValueArgs.number)
        numberOfInterpolationPoints = NameValueArgs.number;
        interpolationIncrement = [];
    elseif isempty(NameValueArgs.number) && isempty(NameValueArgs.increment)
        interpolationIncrement = 1e-3;
        numberOfInterpolationPoints = [];
    end
    sTolerance = NameValueArgs.sTolerance;
    maxIterations = NameValueArgs.maxIterations;
    interpolationMethod = NameValueArgs.interpolationMethod;
    %
    %% create expanded vector
    %
    nVar = size(varData,1);
    xData = [varData; paraData];
    %
    %% get last solution point
    if numel(xData(1,:))>1
        XLast = xData(:,end);
    else
        XLast = [];
    end
    %% initical calculation of stepsize
    %
    if isempty(sData)
        sData = [0,cumsum(sqrt(sum((xData(:,2:end) - xData(:,1:end-1)).^2,1)))];
    end
    %
    %% iteration process of interpolation
    %
    diffS = inf;
    kIt = 0;
    %
    while diffS > sTolerance && kIt < maxIterations
        if ~isempty(numberOfInterpolationPoints)
            sInterp = linspace(sData(1), sData(end), numberOfInterpolationPoints);
        else
            sInterp = sData(1):interpolationIncrement:sData(end);
        end
        %
        % Calculate interpolated points of curve
        xInterp = (interp1(sData, xData.', sInterp, interpolationMethod)).';
        %
        % With those recalculate arclength
        sInterp = [0,cumsum(sqrt(sum((xInterp(:,2:end) - xInterp(:,1:end-1)).^2,1)))];
        %
        % caculate difference in covered arclength
        diffS = abs(sData(end) - sInterp(end));
        %
        % update
        sData = sInterp;
        xData = xInterp;
        %
        kIt = kIt + 1;
    end
    % add Last point
    xInterp = [xInterp,XLast];
    sInterp = [sInterp,norm(diff(xInterp(:,end+(-1:0)),1,2))];
    %% create output
    varInterp = xInterp(1:nVar,:);
    paraInterp = xInterp((nVar+1):end,:);
end

function mustBeIncreasing(array)
    if any(diff(array) <= 0)
        eidType = 'smoothenPath:notIncreasing';
        msgType = 'Input must be an array of increasing values.';
        error(eidType,msgType)
    end
end