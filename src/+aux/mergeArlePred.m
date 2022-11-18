%% path continuation - aux.mergeArlePred
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   24.11.2020 - Alwin Förster
%
function [fs,Js] = mergeArlePred(funPredictor,resCorr,s,xi,ds,Jac)
    [xp,dxpds] = funPredictor(s);
    [ra,dradxp] = resCorr(xp,xi,ds,Jac);
    fs = ra;
    Js = dradxp*dxpds;
end