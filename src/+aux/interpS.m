%% path continuation - aux.interpS
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   27.10.2020 - Tido Kubatschek
%
function [sAllI] = interpS(sAll, varAll, lAll, Opt)
    %%
    Opt.interpolationMethod.spline = true;
    Opt.numberInterpS = 100;
    numberOfPoints = length(lAll);
    
    if numberOfPoints < 2
        error('for interpolation at least 2 data points are needed.')
    end
    
    xAll = [varAll; lAll];
    
    if Opt.interpolationMethod.spline
        % use MATLAB intern function spline
        pp = @(s,x) spline(s, x);

    elseif Opt.interpolationMethod.pchip
        % use MATLAB intern function pchip
        pp = @(s,x) pchip(s, x);
        
    elseif Opt.interpolationMethod.makima
        % use MATLAB intern function makima
        pp = @(s,x) makima(s, x);
        
    else
        error('There is no such interpolation method.')
    end
    
    % calc new s with finer inc
    
    numberOfPoints = min([numberOfPoints, 4]);
    sOld = 0;
    sStart = sAll(end-numberOfPoints);
    %%
    sAllTmp = sAll;
    while true
        sI = linspace(sAllTmp(end-numberOfPoints), sAllTmp(end), Opt.numberInterpS);
        xAllI = ppval(pp(sAllTmp, xAll), sI);

        sAllN = sStart;

        for k = 2:length(xAllI)
            sAllN = [sAllN,sAllN(end) + norm(xAllI(:,k) - xAllI(:,k-1))];
        end
        
        if abs(sAllN(end) - sOld) <= 0.1
            break;
        end
        
        sOld = sAllN(end);
        xAll = xAll(:,1:end);
        sAllTmp = sAllTmp(1:end);
        
        xAll = [xAll, xAllI(:,2:end)];
        sAllTmp = [sAllTmp, sAllN(2:end)];
    end
    sAllTmp = sAllTmp(1:end);
    sAllTmp = [sAllTmp, sAllN(2:end)];
    sAllI = sAllTmp;
    
    %%
end