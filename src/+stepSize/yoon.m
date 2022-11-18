%% path continuation - stepSize.yoon
%  Adjusts stepsize by the ratio of the angles of the lines connecting the
%  last four consecutive solution points. Also adapts due to needed number 
%  of iterations (see <a href="matlab:doc('stepSize.iterationsPolynomial')">stepSize.iterationsPolynomial</a>).
%  Both adaption factors are weighted by the weights specified in 
%  'weightsYoon' and then multiplied.
%
%
%   Inputs:
%       Solver.output -- contains information of solver, such as the 
%                        needed number of iterations.
%       Path          -- contains the solution points of the path
%       Opt           -- contains user inputs, such as optimal number of
%                        iterations, accessible by 'nIterOpt' and the
%                        weights, accessible by 'weightsYoon'.
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('stepSize.control')">other stepsize adaption methods</a>.
%
%   DOI: 10.1177/0954406215586588  (adapted version)
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.01.2022 - Tido Kubatschek
%
function [xi] = yoon(Solver,Path,Opt)
    %% Method of Yoon and Kim
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
    % adapt stepsize
    %
    xiU = 1.1; xiL = 1e10;
    a = (xiU - 1) / (1 - xiL);
    xi = (xiU - ((xiU - xiL) * a) / ((angleM1/angle) + a));
    %
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
    % correct stepsize by iterations
    %
    xi = xi * deviationOfIterations^Opt.weightsYoon;
    %
end