%% path continuation - aux.update_Opt
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.01.2022 - Alwin Förster
%
function [Opt] = update_Opt(Opt,Opt_is_set,Info)
    %% set direction if l_0~=l_start || l_target~=l_end
    %
    if ~Opt_is_set.direction
        if Opt_is_set.l_0
            if Opt_is_set.l_target
                Opt.direction = sign(Opt.l_target-Opt.l_0)*[zeros(size(Info.var0));1];
            else
                Opt.direction = sign(Info.l_end-Opt.l_0)*[zeros(size(Info.var0));1];
            end
        else
            if Opt_is_set.l_target
                Opt.direction = sign(Opt.l_target-Info.l_start)*[zeros(size(Info.var0));1];
            end
        end
    end
    %
    %% set ds_tol dependent on corrector method:
    %
    if ~Opt_is_set.ds_tol
        if Opt.corrector.sphere
            Opt.ds_tol = [0.99,1.01];
        elseif Opt.corrector.orthogonal
            Opt.ds_tol = [0.99,5];
        elseif Opt.corrector.ellipsoid
            Opt.ds_tol = [0.24,1.01];
        elseif Opt.corrector.ellipsoid2
            Opt.ds_tol = [0.000001,1.01];
        elseif Opt.corrector.unique
            Opt.ds_tol = [0.99,5];
        elseif Opt.corrector.paraboloid
            Opt.ds_tol = [0.09,1.01];
        else
            Opt.ds_tol = [0.5,1.5];
        end
    end
    %
end