%% path continuation - stepSize.fayezioghani
%  Adjusts stepsize by the ratio of the angle of the lines connecting the
%  last three consecutive solution points and the tangent in the current 
%  point. Also adapts due to needed number of iterations (see 
%  <a href="matlab:doc('stepSize.iterationsPolynomial')">stepSize.iterationsPolynomial</a>).
%  Both adaption factors are weighted by the weights specified in 
%  'weightsFayezioghani' and then multiplied.
%
%
%   Inputs:
%       ds            -- latest used stepsize
%       Solver.output -- contains information of solver, such as the 
%                        needed number of iterations.
%       Path          -- contains the solution points of the path
%       Jac           -- current Jacobian to calculate tangent
%       Opt           -- contains user inputs, such as optimal number of
%                        iterations, accessible by 'nIterOpt' and the
%                        weights, accessible by 'weightsFayezioghani'.
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('stepSize.control')">other stepsize adaption methods</a>.
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
    varNeeded = Path.varAll(:,end-1:end);
    lNeeded = Path.lAll(end-1:end);
    zNeeded = [varNeeded; lNeeded];
    %
    % calculate connecting vector and tangent
    %
    v = zNeeded(:,end) - zNeeded(:,end-1);
    if diff(size(Jac)) == 0 && size(Jac,1) == (size(varNeeded,1) + 1)
        [~,tangent] = predictor.ode(Path,ds,Jac,[],Opt);
    else
        tangent = v;
    end
    %
    % calculate angle
    %
    angle = aux.vectorAngle(v,tangent);
    %
    % get weigths
    %
    weights = Opt.weightsFayezioghani;
    %
    % correct number of iterations
    %
    if Opt.dsMax==inf
        iter = max(Solver.output.iterations(end),1);
    else
        iter = Solver.output.iterations(end);
    end
    %
    % calculate new step size
    %
    xi = (Opt.nIterOpt / iter)^weights(1)*...
        ((cos(angle) + 1) / (cos(Opt.stepSizeAngle) + 1))^weights(2);
end