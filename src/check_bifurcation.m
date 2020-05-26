%% path continuation - check_bifurcation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2020 - Alwin Förster
%
function [var_all,l_all,bif,sign_det_jacobian] = check_bifurcation(fun,solver_jacobian_red,var_all,l_all,bif,sign_det_jacobian,Opt)
    if Opt.bifurcation.mark
        %% mark bifurcations:
        full_rank = length(solver_jacobian_red(:,1));
        rank_tol = 10^-2; % TODO!!!
        sign_det_current_jacobian = sign(det(solver_jacobian_red));
        if sign_det_current_jacobian*sign_det_jacobian<0
            bifurcation_type = (rank(solver_jacobian_red,rank_tol)==full_rank); % 1: true bif.; 0: zero point
            bif = [bif,[length(l_all);bifurcation_type]];
            sign_det_jacobian = sign_det_current_jacobian;
        end
    elseif Opt.bifurcation.determine || Opt.bifurcation.trace
        %% determine bifurcation-points
        warning('not implemented, using mark instead.');
        Opt.bifurcation.mark = true;
        Opt.bifurcation.determine = false;
        Opt.bifurcation.trace = false;
        %% mark bifurcations:
        full_rank = length(solver_jacobian_red(:,1));
        rank_tol = 10^-2; % TODO!!!
        sign_det_current_jacobian = sign(det(solver_jacobian_red));
        if sign_det_current_jacobian*sign_det_jacobian<0
            bifurcation_type = (rank(solver_jacobian_red,rank_tol)==full_rank); % 1: true bif.; 0: zero point
            bif = [bif,[length(l_all);bifurcation_type]];
            sign_det_jacobian = sign_det_current_jacobian;
        end
    else
        %% error
        error('No such bifurcation-mode');
    end
end