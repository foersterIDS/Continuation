%% path continuation - shorten_dummy
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.02.2021 - Tido Kubatschek
%
function [x_solution, fun_solution, solver_jacobian, dscale, x_predictor] = shorten_dummy_results(x_solution, fun_solution, solver_jacobian, dscale, x_predictor)
    x_solution = x_solution(1:(end-1));
    fun_solution = fun_solution(1:(end-1));
    solver_jacobian = solver_jacobian(1:(end-1), 1:(end-1));
    dscale = dscale(1:end-1);
    x_predictor = x_predictor(1:end-1);
end