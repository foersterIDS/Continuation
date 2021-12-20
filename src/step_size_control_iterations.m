%% path continuation - step_size_control_iterations
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.06.2020 - Niklas Marhenke
%   08.10.2020 - Alwin Förster
%
%   DOI: 10.1007/978-1-4419-1740-9
%
function [dsn] = step_size_control_iterations(ds,ds0,Counter,solver_output,Do,Path,Opt)
    % calculate step size
    if ison(Opt.step_size_iterations_type.polynomial)
        dsn = ds*(Opt.n_iter_opt/solver_output.iterations)^Opt.step_size_iterations_beta;
    else
        dsn = ds*2^((solver_output.iterations - Opt.n_iter_opt)/4);
    end
end