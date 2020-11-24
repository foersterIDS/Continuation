%% path continuation - merge_arle_pred
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   24.11.2020 - Alwin Förster
%
function [fs,Js] = merge_arle_pred(fun_predictor,res_arle,s,xi,ds)
    [xp,dxpds] = fun_predictor(s);
    [ra,dradxp] = res_arle(xp,xi,ds);
    fs = ra;
    Js = dradxp*dxpds;
end