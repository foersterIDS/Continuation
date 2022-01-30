%% path continuation - step_size.pid_valli
%  Adjusts stepsize by the relative change of the variables which is
%  implemented in a PID control. The control constants of the P, I and D
%  controller can be specified by 'step_size_pid_params'. The relative
%  change is calculated with respect to a tolerance specified in 
%  'step_size_pid_tol'. Needed derivatives are calculated with respect to
%  the varying parameter lamdba.
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
%  Also see <a href="matlab:doc('step_size.pid_custom')">step_size.pid_custom</a> or
%  see the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('step_size.control')">other stepsize adaption methods</a>.
%
%   DOI: 10.1002/ﬂd.1998
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   21.10.2020 - Tido Kubatschek
% 
function [xi] = pid_valli(Path,Opt)
    % define function handle for finite difference calculation
    %
    dvarsdl = @(k) (Path.var_all(:,end+k) - Path.var_all(:,end+k-1)) / (Path.l_all(end+k) - Path.l_all(end+k-1));
    %
    % define function handle for e_tilde
    %
    e_tilde= @(k) norm(dvarsdl(k) - dvarsdl(k-1)) / norm(dvarsdl(k));
    %
    % parameters of pid control 
    %
    PID = Opt.step_size_pid_params;
    pid_tol = Opt.step_size_pid_tol;
    %
    % calculate needed relative changes
    %
    e = [e_tilde(0), e_tilde(-1), e_tilde(-2)];
    e = e / pid_tol;
    %
    % calculate step size adaption factor
    %
    xi = (e(2) / e(1)) ^ PID(1) *...
        (1 / e(1)) ^ PID(2) *...
        (e(2)^2 / (e(1) * e(3)) ) ^ PID(3);
    %
end