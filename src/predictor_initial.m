%% path continuation - predictor_initial
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   24.11.2020 - Alwin Förster
%
function [fun_predictor,Jac_predictor] = predictor_initial(Path,s,Opt)
    if numel(Opt.direction)==1
        fun_predictor = [Path.var_all(:,end);Path.l_all(end)+sign(Opt.direction)*s];
        Jac_predictor = [zeros(size(Path.var_all(:,end)));sign(Opt.direction)];
    else
        fun_predictor = [Path.var_all(:,end);Path.l_all(end)]+Opt.direction*s;
        Jac_predictor = Opt.direction;
    end
end