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
%       oih           -- OptInfoHandle object
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
function [xi] = fayezioghani(ds,oih)
    %% Method of Fayezioghani et al.
    %
    % collect needed data of Path
    %
    varNeeded = oih.path.varAll(:,end-1:end);
    lNeeded = oih.path.lAll(end-1:end);
    zNeeded = [varNeeded; lNeeded];
    %
    % calculate connecting vector and tangent
    %
    jac = oih.path.getJacobianByName('last');
    v = zNeeded(:,end) - zNeeded(:,end-1);
    if diff(size(jac)) == 0 && size(jac,1) == (size(varNeeded,1) + 1)
        [~,tangent] = predictor.ode(oih,ds,jac,[]);
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
    weights = oih.opt.weightsFayezioghani;
    %
    % correct number of iterations
    %
    if oih.opt.dsMax==inf
        iter = max(oih.solver.output.iterations(end),1);
    else
        iter = oih.solver.output.iterations(end);
    end
    %
    % calculate new step size
    %
    xi = (oih.opt.nIterOpt / iter)^weights(1)*...
        ((cos(angle) + 1) / (cos(oih.opt.stepSizeAngle) + 1))^weights(2);
end