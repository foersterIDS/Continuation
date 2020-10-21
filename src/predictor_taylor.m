%% path continuation - predictor_taylor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   05.10.2020 - Alwin Förster
%
function [xp] = predictor_taylor(var_all,l_all,s_all,no,nf,ds)
    no = min([length(l_all)-1,no]);
    ns = min([length(l_all),no+1+nf]);
    xs = [var_all;l_all];
    nd = length(xs(:,end));
    %% calc taylor-predictor:
    p = polyfitn(s_all(end+((-ns+1):0)),xs(:,end+((-ns+1):0)),no);
    xp = polyvaln(p,s_all(end)+ds,nd);
end