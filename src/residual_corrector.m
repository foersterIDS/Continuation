%% path continuation - residual_corrector
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   17.09.2020 - Tido Kubatschek
%
function [residual] = residual_corrector(Opt)
    if Opt.corrector.linear
        residual = @(x,xs,ds) residual_corrector_linear(x,xs,ds,Opt);
    elseif Opt.corrector.sphere
        residual = @(x,xs,ds) residual_corrector_sphere(x,xs,ds);
    elseif Opt.corrector.ellipsoid
        RnR = RnRotation([1;0]);
        residual = @(x,xs,ds) residual_corrector_ellipsoid(x,xs,ds,RnR);
    elseif Opt.corrector.ellipsoid2
        residual = @(x,xs,ds) residual_corrector_ellipsoid2(x,xs,ds);
    elseif Opt.corrector.unique
        residual = @(x,xs,ds) residual_corrector_unique(x,xs,ds,Opt);
    else
        error('No such corrector-method');
    end
end