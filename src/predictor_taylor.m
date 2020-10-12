%% path continuation - predictor_taylor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   05.10.2020 - Alwin Förster
%
function [xp] = predictor_taylor(vars,ls,no,nf,ds)
    no = min([length(ls)-1,no]);
    ns = min([length(ls),no+1+nf]);
    xs = [vars;ls];
    nd = length(xs(:,end));
    %% calc arc-length approximation:
    s = zeros(1,ns);
    for i=2:ns
        s(i) = s(i-1)+norm(xs(:,end-ns+i)-xs(:,end-ns+i-1));
    end
    %% calc taylor-predictor:
    p = polyfitn(s,xs(:,end+((-ns+1):0)),no);
    xp = polyvaln(p,s(end)+ds,nd);
end