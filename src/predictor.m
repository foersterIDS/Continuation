%% path continuation - predictor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [vp,lp] = predictor(vars,ls,ds,Opt)
    if length(ls)==1
        xip1 = [vars;ls+ds];
    else
        xi = [vars(:,end);ls(end)];
        xim1 = [vars(:,end-1);ls(end-1)];
        xip1 = xi+(xi-xim1)*ds/sqrt((xi-xim1)'*(xi-xim1));
    end
    vp = xip1(1:end-1);
    lp = xip1(end);
end