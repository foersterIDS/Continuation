%% path continuation - stepSize.iterationsExponential
%  Adjusts stepsize by needed number of iterations due to the difference of 
%  optimal number and the needed number of iterations specified by the
%  user. The difference then used as an exponent and can be weighted by
%  'stepSizeExponentialWeight'.
%
%
%   Inputs:
%       oih           -- OptInfoHandle object
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  Also see <a href="matlab:doc('stepSize.iterationsPolynomial')">stepSize.iterationsPolynomial</a> or
%  see the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('stepSize.control')">other stepsize adaption methods</a>.
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.10.2020 - Alwin Förster
%   17.01.2022 - Tido Kubatschek
%
function [xi] = iterationsExponential(oih)
    % correct number of iterations
    if oih.opt.dsMax==inf
        iter = max(oih.solver.output.iterations(end),1);
    else
        iter = oih.solver.output.iterations(end);
    end
    % calculate step size adaption factor
    xi = 2^((oih.opt.nIterOpt - iter)/oih.opt.stepSizeExponentialWeight);
end