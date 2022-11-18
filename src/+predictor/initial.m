%% path continuation - predictor.initial
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   24.11.2020 - Alwin Förster
%
function [funPredictor,JacPredictor] = initial(Path,s,Opt)
    if numel(Opt.direction)==1
        funPredictor = [Path.varAll(:,end);Path.lAll(end)+sign(Opt.direction)*s];
        JacPredictor = [zeros(size(Path.varAll(:,end)));sign(Opt.direction)];
    else
        funPredictor = [Path.varAll(:,end);Path.lAll(end)]+Opt.direction*s;
        JacPredictor = Opt.direction;
    end
end