%% path continuation - add_dummy
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.02.2021 - Tido Kubatschek
%
function [residual, x_predictor, dscale] = add_dummy(residual, ds, var_all, l_all, Opt, x_predictor, dscale)
    % adjust residual
    residual = @(x) merge_residuals(@(x,g) residual(x),@(x,xs,ds) residual_dummy(x,xs,ds), x, [var_all;l_all],ds,Opt);
    % adjust predictor and dscale
    x_predictor = [x_predictor; 1.001];
    dscale = [dscale;1];
end