%% path continuation - step_size.fayezioghani
%  Adjusts stepsize by the ratio of the angle of the lines connecting the
%  last three consecutive solution points and the tangent in the current 
%  point. Also adapts due to needed number of iterations (see 
%  <a href="matlab:doc('step_size.iterations_polynomial')">step_size.iterations_polynomial</a>).
%  Both adaption factors are weighted by the weights specified in 
%  'weights_fayezioghani' and then multiplied.
%
%
%   Inputs:
%       ds            -- latest used stepsize
%       Solver.output -- contains information of solver, such as the 
%                        needed number of iterations.
%       Path          -- contains the solution points of the path
%       Jac           -- current Jacobian to calculate tangent
%       Opt           -- contains user inputs, such as optimal number of
%                        iterations, accessible by 'n_iter_opt' and the
%                        weights, accessible by 'weights_fayezioghani'.
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('step_size.control')">other stepsize adaption methods</a>.
%
%   DOI: https://doi.org/10.1016/j.compstruc.2019.07.009 (adapted version)
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.01.2022 - Tido Kubatschek
%
function [xi] = fayezioghani(ds,Solver,Path,Jac,Opt)
    %% Method of Fayezioghani et al.
    %
    % collect needed data of Path
    %
    var_needed = Path.var_all(:,end-1:end);
    l_needed = Path.l_all(end-1:end);
    z_needed = [var_needed; l_needed];
    %
    % calculate connecting vector and tangent
    %
    v = z_needed(:,end) - z_needed(:,end-1);
    if diff(size(Jac)) == 0 && size(Jac,1) == (size(var_needed,1) + 1)
        [~,tangent] = predictor.ode(Path,ds,Jac,[],Opt);
    else
        tangent = v;
    end
    %
    % calculate angle
    %
    angle = aux.vector_angle(v,tangent);
    %
    % get weigths
    %
    weights = Opt.weights_fayezioghani;
    %
    % correct number of iterations
    %
    if Opt.ds_max==inf
        iter = max(Solver.output.iterations(end),1);
    else
        iter = Solver.output.iterations(end);
    end
    %
    % calculate new step size
    %
    xi = (Opt.n_iter_opt / iter)^weights(1)*...
        ((cos(angle) + 1) / (cos(Opt.step_size_angle) + 1))^weights(2);
end