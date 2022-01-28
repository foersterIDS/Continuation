%% path continuation - step_size.pid_custom
%  Adjusts stepsize by the relative change of the variables which is
%  implemented in a PID control. The control constants of the P, I and D
%  controller can be specified by 'step_size_pid_params'. The relative
%  change is calculated with respect to a tolerance specified in 
%  'step_size_pid_tol'. This method is inspired by an adaption method
%  mentioned in the literature below. In contrast this method calculates
%  needed derivatives with respect to arclength.
%
%
%   Inputs:
%       Path    -- contains information the path, such as the the values of
%                  lamdba and the variables
%       Opt     -- contains user inputs, such as PID constants, accessible
%                  by 'step_size_pid_params' and the tolerance, accessible
%                  by 'step_size_pid_tol'.
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  Also see <a href="matlab:doc('step_size.pid_valli')">step_size.pid_valli</a> or
%  see the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('step_size.control')">other stepsize adaption methods</a>.
%
%   DOI: 10.1002/ﬂd.1998 (adapted version)
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.01.2022 - Tido Kubatschek
%
function [xi] = pid_custom(Path,Opt)
    z_all = [Path.var_all; Path.l_all];
    dzds = @(k) (z_all(:,end+k) - z_all(:,end+k-1)) / (Path.s_all(end+k) - Path.s_all(end+k-1));
    e_star = @(k) norm(dzds(k) - dzds(k-1)) / norm(dzds(k));
    %
    % parameters of pid control 
    %
    PID = Opt.step_size_pid_params;
    pid_tol = Opt.step_size_pid_tol;
    %
    % relative change
    %
    e = [e_star(0), e_star(-1), e_star(-2)];
    e = e / pid_tol;
    %
    % calculate step size
    %
    xi = (e(2) / e(1))^PID(1) * (1 / e(1))^PID(2) * (e(2)^2 / (e(1) * e(3)) )^PID(3);
    %
end