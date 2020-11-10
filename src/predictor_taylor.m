%% path continuation - predictor_taylor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   05.10.2020 - Alwin Förster
%
function [xp] = predictor_taylor(var_all,l_all,s_all,no,nf,ds)
    no = min([length(l_all)-1,no]);
    ns = min([length(l_all),no+1+nf]);
    x_all = [var_all;l_all];
    nd = length(x_all(:,end));
    %% calc scaling:
    dsc_min = 10^-15;
    x_basis = x_all(:,end-ns+1);
    s_basis = s_all(end-ns+1);
    if length(l_all)>2
        dsc_x = max([mean(diff(abs(x_all(:,end+((-ns+1):0))-x_basis),1,2),2),ones(nd,1)*dsc_min]')';
        dsc_s = max([mean(diff(abs(s_all(end+((-ns+1):0))-s_basis),1,2),2),dsc_min]);
    else
        dsc_x = max([abs(x_all(:,2)-x_basis),ones(nd,1)*dsc_min]')';
        dsc_s = max([abs(s_all(2)-s_basis),dsc_min]);
    end
    %% calc taylor-predictor:
    p_sc = polyfitn((s_all(end+((-ns+1):0))-s_basis)./dsc_s,(x_all(:,end+((-ns+1):0))-x_basis)./kron(dsc_x,ones(1,ns)),no);
    xp_sc = polyvaln(p_sc,(s_all(end)+ds-s_basis)/dsc_s,nd);
    %% descale predictor:
    xp = xp_sc.*dsc_x+x_basis;
end