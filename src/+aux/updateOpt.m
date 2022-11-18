%% path continuation - aux.updateOpt
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.01.2022 - Alwin Förster
%
function [Opt] = updateOpt(Opt,OptIsSet,Info)
    %% set direction if l0~=lStart || lTarget~=lEnd
    %
    if ~OptIsSet.direction
        if OptIsSet.l0
            if OptIsSet.lTarget
                Opt.direction = sign(Opt.lTarget-Opt.l0)*[zeros(size(Info.var0));1];
            else
                Opt.direction = sign(Info.lEnd-Opt.l0)*[zeros(size(Info.var0));1];
            end
        else
            if OptIsSet.lTarget
                Opt.direction = sign(Opt.lTarget-Info.lStart)*[zeros(size(Info.var0));1];
            end
        end
    end
    %
    %% set dsTol dependent on corrector method:
    %
    if ~OptIsSet.dsTol
        if Opt.corrector.sphere
            Opt.dsTol = [0.99,1.01];
        elseif Opt.corrector.orthogonal
            Opt.dsTol = [0.99,5];
        elseif Opt.corrector.ellipsoid
            Opt.dsTol = [0.24,1.01];
        elseif Opt.corrector.ellipsoid2
            Opt.dsTol = [0.000001,1.01];
        elseif Opt.corrector.unique
            Opt.dsTol = [0.99,5];
        elseif Opt.corrector.paraboloid
            Opt.dsTol = [0.09,1.01];
        else
            Opt.dsTol = [0.5,1.5];
        end
    end
    %
end