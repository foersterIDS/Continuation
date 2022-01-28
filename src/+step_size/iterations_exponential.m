%% path continuation - step_size.iterations_exponential
%  Adjusts stepsize by needed number of iterations due to the difference of 
%  optimal number and the needed number of iterations specified by the
%  user. The difference then used as an exponent and can be weighted by
%  'step_size_exponential_weight'.
%
%
%   Inputs:
%       solver_output -- contains information of solver, such as the 
%                        needed number of iterations.
%       Opt           -- contains user inputs, such as optimal number of
%                        iterations, accessible by 'n_iter_opt' and the
%                        exponent, accessible by 
%                        'step_size_exponential_weight'
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  Also see <a href="matlab:doc('step_size.iterations_polynomial')">step_size.iterations_polynomial</a> or
%  see the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('step_size.control')">other stepsize adaption methods</a>.
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.10.2020 - Alwin Förster
%   17.01.2022 - Tido Kubatschek
%
function [xi] = iterations_exponential(solver_output,Opt)
    % correct number of iterations
    if Opt.ds_max==inf
        iter = max(solver_output.iterations(end),1);
    else
        iter = solver_output.iterations(end);
    end
    % calculate step size adaption factor
    xi = 2^((Opt.n_iter_opt - iter)/Opt.step_size_exponential_weight);
end