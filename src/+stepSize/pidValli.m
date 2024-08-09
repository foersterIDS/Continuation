%% path continuation - stepSize.pidValli
%  Adjusts stepsize by the relative change of the variables which is
%  implemented in a PID control. The control constants of the P, I and D
%  controller can be specified by 'stepSizePidParams'. The relative
%  change is calculated with respect to a tolerance specified in 
%  'stepSizePidTol'. Needed derivatives are calculated with respect to
%  the varying parameter lamdba.
%
%
%   Inputs:
%       oih           -- OptInfoHandle object
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  Also see <a href="matlab:doc('stepSize.pidCustom')">stepSize.pidCustom</a> or
%  see the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('stepSize.control')">other stepsize adaption methods</a>.
%
%   DOI: 10.1002/ﬂd.1998
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   21.10.2020 - Tido Kubatschek
% 
function [xi] = pidValli(oih)
    % define function handle for finite difference calculation
    %
    dvarsdl = @(k) (oih.path.varAll(:,end+k) - oih.path.varAll(:,end+k-1)) / (oih.path.lAll(end+k) - oih.path.lAll(end+k-1));
    %
    % define function handle for eTilde
    %
    eTilde= @(k) norm(dvarsdl(k) - dvarsdl(k-1)) / norm(dvarsdl(k));
    %
    % parameters of pid control 
    %
    PID = oih.opt.stepSizePidParams;
    pidTol = oih.opt.stepSizePidTol;
    %
    % calculate needed relative changes
    %
    e = [eTilde(0), eTilde(-1), eTilde(-2)];
    e = e / pidTol;
    %
    % calculate step size adaption factor
    %
    xi = (e(2) / e(1)) ^ PID(1) *...
        (1 / e(1)) ^ PID(2) *...
        (e(2)^2 / (e(1) * e(3)) ) ^ PID(3);
    %
end