%% path continuation - predictor_initial
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   24.11.2020 - Alwin Förster
%
function [fun_predictor,Jac_predictor] = predictor_initial(var_all,l_all,s,Opt)
    if numel(Opt.direction)==1
        fun_predictor = [var_all(:,end);l_all(end)+sign(Opt.direction)*s];
        Jac_predictor = [zeros(size(var_all(:,end)));sign(Opt.direction)];
    else
        fun_predictor = [var_all(:,end);l_all(end)]+Opt.direction*s;
        Jac_predictor = Opt.direction;
    end
end