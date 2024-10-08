%% path continuation - stepSize.pidCustom
%  Adjusts stepsize by the relative change of the variables which is
%  implemented in a PID control. The control constants of the P, I and D
%  controller can be specified by 'stepSizePidParams'. The relative
%  change is calculated with respect to a tolerance specified in 
%  'stepSizePidTol'. This method is inspired by an adaption method
%  mentioned in the literature below. In contrast this method calculates
%  needed derivatives with respect to arclength.
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
%  Also see <a href="matlab:doc('stepSize.pidValli')">stepSize.pidValli</a> or
%  see the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('stepSize.control')">other stepsize adaption methods</a>.
%
%   DOI: 10.1002/ﬂd.1998 (adapted version)
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.01.2022 - Tido Kubatschek
%
function [xi] = pidCustom(oih)
    % solution points
    %
    zAll = oih.path.xAll;
    %
    % define function handle for finite difference calculation
    %
    dzds = @(k) (zAll(:,end+k) - zAll(:,end+k-1)) / (oih.path.sAll(end+k) - oih.path.sAll(end+k-1));
    %
    % define function handle for eTilde
    %
    eTilde = @(k) norm(dzds(k) - dzds(k-1)) / norm(dzds(k));
    %
    % parameters of pid control 
    %
    PID = oih.opt.stepSizePidParams;
    pidTol = oih.opt.stepSizePidTol;
    %
    % relative change
    %
    e = [eTilde(0), eTilde(-1), eTilde(-2)];
    e = e / pidTol;
    %
    % calculate stepsize adaption factor
    %
    xi = (e(2) / e(1))^PID(1) * (1 / e(1))^PID(2) * (e(2)^2 / (e(1) * e(3)) )^PID(3);
    %
end