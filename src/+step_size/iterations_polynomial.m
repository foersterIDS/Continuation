%% path continuation - step_size.iterations_polynomial
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.06.2020 - Niklas Marhenke
%   08.10.2020 - Alwin Förster
%   17.01.2022 - Tido Kubatschek
%
%   DOI: 10.1007/978-1-4419-1740-9
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


