%% path continuation - aux.initialBreakFunction
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   21.02.2021 - Alwin Förster
%
function [doBreak,breakFunOut] = initialBreakFunction(fun,jac,v,l,breakFunOut)
    doBreak = false;
    breakFunOut = [breakFunOut,numel(breakFunOut)+1];
end