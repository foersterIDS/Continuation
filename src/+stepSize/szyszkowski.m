%% path continuation - stepSize.szyszkowski
%  Adjusts stepsize by the ratio of the angles of the lines connecting the
%  last four consecutive solution points. Also adapts due to needed number 
%  of iterations (see <a href="matlab:doc('stepSize.iterationsPolynomial')">stepSize.iterationsPolynomial</a>).
%  Both adaption factors are weighted by the weights specified in 
%  'weightsSzyszkowski' and then multiplied.
%
%
%   Inputs:
%       Solver.output -- contains information of solver, such as the 
%                        needed number of iterations.
%       Path          -- contains the solution points of the path
%       Opt           -- contains user inputs, such as optimal number of
%                        iterations, accessible by 'nIterOpt' and the
%                        weights, accessible by 'weightsSzyszkowski'.
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('stepSize.control')">other stepsize adaption methods</a>.
%
%   DOI: https://doi.org/10.1007/s004660050513  (adapted version)
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.01.2022 - Tido Kubatschek
%
function [xi] = szyszkowski(Solver,Path,Opt)
    %% Method of Szyszkowski and Husband
    %
    % collect needed data of Path
    %
    varNeeded = Path.varAll(:,end-3:end);
    lNeeded = Path.lAll(end-3:end);
    zNeeded = [varNeeded; lNeeded];
    %
    % calculate connecting vectors
    %
    v1 = zNeeded(:,end) - zNeeded(:,end-1);
    v2 = zNeeded(:,end-1) - zNeeded(:,end-2);
    v3 = zNeeded(:,end-2) - zNeeded(:,end-3);
    %
    % calculate angles
    %
    angle = aux.vectorAngle(v1,v2);
    angleM1 = aux.vectorAngle(v2,v3);
    %
    % calculate next angle
    %
    angleP1 = 2*angle - angleM1 * norm(v1)/norm(v2);
    %
    % check if angle is too large
    %
    if angle > Opt.stepSizeAngle || angleP1 > Opt.stepSizeAngle
        xi = 0.5;
    else
        % correct number of iterations
        %
        if Opt.dsMax==inf
            iter = max(Solver.output.iterations(end),1);
        else
            iter = Solver.output.iterations(end);
        end
        %
        % calculate deviation of iterations
        %
        deviationOfIterations = Opt.nIterOpt/iter;
        %
        % get weigths
        weights = Opt.weightsSzyszkowski;
        %
        % calculate new step size
        %
        xi = deviationOfIterations^weights(1) *...
            abs(angle/angleP1)^weights(2);
        %
    end
end