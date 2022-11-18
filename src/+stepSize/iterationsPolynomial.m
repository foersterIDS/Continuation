%% path continuation - stepSize.iterationsPolynomial
%  Adjusts stepsize by needed number of iterations due to the quotient of 
%  needed number and the optimal number of iterations specified by the
%  user. The quotient is then raised to the power of 
%  'stepSizeIterationsBeta'.
%
%
%   Inputs:
%       Solver.output -- contains information of solver, such as the 
%                        needed number of iterations.
%       Opt           -- contains user inputs, such as optimal number of
%                        iterations, accessible by 'nIterOpt' and the
%                        exponent, accessible by 'stepSizeIterationsBeta'
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
function [xi] = iterationsPolynomial(Solver,Opt)
    % correct number of iterations
    if Opt.dsMax==inf
        iter = max(Solver.output.iterations(end),1);
    else
        iter = Solver.output.iterations(end);
    end
    % calculate step size adaption factor
    xi = (Opt.nIterOpt/iter)^Opt.stepSizeIterationsBeta;
end


