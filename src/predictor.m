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
        [nt,nf] = predictor_adaptive(vars,ls,Opt);
        xip1 = predictor_taylor(vars,ls,nt,nf,ds);
    end
    vp = xip1(1:end-1);
    lp = xip1(end);
end