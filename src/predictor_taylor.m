%% path continuation - predictor_taylor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   05.10.2020 - Alwin Förster
%
function [fun_predictor,Jac_predictor] = predictor_taylor(var_all,l_all,s_all,no,nf)
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
    fun_predictor_sc = @(s) polyvaln(p_sc,(s_all(end)+s-s_basis)/dsc_s,nd);
    if nargout>1
        dp_sc = polydern(p_sc,nd);
        Jac_predictor_sc = @(s) polyvaln(dp_sc,(s_all(end)+s-s_basis)/dsc_s,nd);
    end
    %% descale predictor:
    fun_predictor = @(s) fun_predictor_sc(s).*dsc_x+x_basis;
    if nargout>1
        Jac_predictor = @(s) Jac_predictor_sc(s).*dsc_x;
    end
end