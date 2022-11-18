%% path continuation - stepSize.iterationsExponential
%  Adjusts stepsize by needed number of iterations due to the difference of 
%  optimal number and the needed number of iterations specified by the
%  user. The difference then used as an exponent and can be weighted by
%  'stepSizeExponentialWeight'.
%
%
%   Inputs:
%       Solver.output -- contains information of solver, such as the 
%                        needed number of iterations.
%       Opt           -- contains user inputs, such as optimal number of
%                        iterations, accessible by 'nIterOpt' and the
%                        exponent, accessible by 
%                        'stepSizeExponentialWeight'
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
function [xi] = iterationsExponential(Solver,Opt)
    % correct number of iterations
    if Opt.dsMax==inf
        iter = max(Solver.output.iterations(end),1);
    else
        iter = Solver.output.iterations(end);
    end
    % calculate step size adaption factor
    xi = 2^((Opt.nIterOpt - iter)/Opt.stepSizeExponentialWeight);
end