%% path continuation - stepSize.contraction
%  Adjusts stepsize by the rate of contraction of the solver which is
%  compared to an optimale rate.
%
%
%   Inputs:
%       Solver.output -- contains information of solver, such as the 
%                        rate of contraction.
%       Opt           -- contains user inputs, such as the optimal contraction
%                        rate specified in 'optimalContractionRate'.
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('stepSize.control')">other stepsize adaption methods</a>.
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.01.2022 - Tido Kubatschek
%
function [xi] = contraction(Solver,Opt)
    if ~isempty(Solver.output.rateOfContraction)
        % get rates of contraction
        %
        % current rate
        %
        currentRate = Solver.output.rateOfContraction(end);
        %
        % if there is more than one rate accessible also get previous rate
        % and calculate relative difference
        %
        if length(Solver.output.rateOfContraction) > 1
            previousRate = Solver.output.rateOfContraction(end-1);
            relDifference = abs(currentRate - previousRate) / abs(previousRate);
        else
            relDifference = inf;
        end
        % 
        % optimal rate
        %
        qOpt = Opt.optimalContractionRate;
        %
        if currentRate > qOpt && relDifference > 0.01
            % only decrease if there was a difference
            xi = 0.5;
        elseif currentRate <= qOpt/4 || (currentRate <= qOpt/2 && relDifference < 0.01)
            % increase if rate is small or there is just a small difference
            xi = 1.2;
        else
            xi = 1;
        end
    else
        xi = 1;
    end
end