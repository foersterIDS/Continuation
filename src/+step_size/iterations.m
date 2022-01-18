%% path continuation - step_size.iterations
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.06.2020 - Niklas Marhenke
%   08.10.2020 - Alwin Förster
%
%   DOI: 10.1007/978-1-4419-1740-9
%
function [dsn] = iterations(ds,ds0,Counter,solver_output,Do,Path,Opt)
    % correct number of iterations
    if Opt.ds_max==inf
        iter = max(solver_output.iterations,1);
    else
        iter = solver_output.iterations;
    end
    % calculate step size
    if aux.ison(Opt.step_size_iterations_type.polynomial)
        dsn = ds*(Opt.n_iter_opt/iter)^Opt.step_size_iterations_beta;
    else
        dsn = ds*2^((iter - Opt.n_iter_opt)/4);
    end
end