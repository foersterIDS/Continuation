%% path continuation - predictor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin F�rster
%
function [vp,lp,ch] = predictor(vars,ls,ds,Opt)
    if length(ls)==1
        xip1 = [vars;ls+sign(Opt.direction)*ds];
        ch=1;
    else
        [xip1,ch] = predictor_taylor(vars,ls,Opt.predictor_taylor,ds);
    end
    vp = xip1(1:end-1);
    lp = xip1(end);
end