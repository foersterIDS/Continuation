%% path continuation - step_size.iterations_polynomial
%  Adjusts stepsize by needed number of iterations due to the quotient of 
%  needed number and the optimal number of iterations specified by the
%  user. The quotient is then raised to the power of 
%  'step_size_iterations_beta'.
%
%
%   Inputs:
%       solver_output -- contains information of solver, such as the 
%                        needed number of iterations.
%       Opt           -- contains user inputs, such as optimal number of
%                        iterations, accessible by 'n_iter_opt' and the
%                        exponent, accessible by 'step_size_iterations_beta'
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  Also see <a href="matlab:doc('step_size.iterations_exponential')">step_size.iterations_exponential</a> or
%  see the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('step_size.control')">other stepsize adaption methods</a>.
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.06.2020 - Niklas Marhenke
%   08.10.2020 - Alwin Förster
%   17.01.2022 - Tido Kubatschek
%
function [xi] = iterations_polynomial(solver_output,Opt)
    % correct number of iterations
    if Opt.ds_max==inf
        iter = max(solver_output.iterations(end),1);
    else
        iter = solver_output.iterations(end);
    end
    % calculate step size adaption factor
    xi = (Opt.n_iter_opt/iter)^Opt.step_size_iterations_beta;
end


