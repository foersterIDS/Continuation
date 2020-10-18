%% path continuation - predictor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [vp,lp] = predictor(vars,ls,ds,Opt)
    if length(ls)==1
        xip1 = [vars;ls+sign(Opt.direction)*ds];
    else
        xip1 = predictor_taylor(vars,ls,Opt.predictor_taylor,Opt.predictor_fit,ds);
    end
    vp = xip1(1:end-1);
    lp = xip1(end);
end