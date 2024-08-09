%% path continuation - predictor.initial
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   24.11.2020 - Alwin Förster
%
function [funPredictor,JacPredictor] = initial(oih,s)
    if numel(oih.opt.direction)==1
        funPredictor = [oih.path.varAll(:,end);oih.path.lAll(end)+sign(oih.opt.direction)*s];
        JacPredictor = [zeros(size(oih.path.varAll(:,end)));sign(oih.opt.direction)];
    else
        funPredictor = [oih.path.varAll(:,end);oih.path.lAll(end)]+oih.opt.direction*s;
        JacPredictor = oih.opt.direction;
    end
end