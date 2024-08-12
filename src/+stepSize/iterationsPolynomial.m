%% path continuation - stepSize.iterationsPolynomial
%  Adjusts stepsize by needed number of iterations due to the quotient of 
%  needed number and the optimal number of iterations specified by the
%  user. The quotient is then raised to the power of 
%  'stepSizeIterationsBeta'.
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
%  Also see <a href="matlab:doc('stepSize.iterationsExponential')">stepSize.iterationsExponential</a> or
%  see the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('stepSize.control')">other stepsize adaption methods</a>.
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.06.2020 - Niklas Marhenke
%   08.10.2020 - Alwin Förster
%   17.01.2022 - Tido Kubatschek
%
function [xi] = iterationsPolynomial(oih)
    % correct number of iterations
    if oih.opt.dsMax==inf
        iter = max(oih.path.iterations(end),1);
    else
        iter = oih.path.iterations(end);
    end
    % calculate step size adaption factor
    xi = (oih.opt.nIterOpt/iter)^oih.opt.stepSizeIterationsBeta;
end


