%% path continuation - homotopy.continuation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [ var_h, exitflag ] = h_continuation( fun, var0, Opt )
    %% Display:
    %
    aux.print_line(Opt,'------> starting homotopy...');
    %
    %% Homotopy-function
    %
    if Opt.homotopy.fix
        G = @(x) homotopy.fix(x,var0);
    elseif Opt.homotopy.newton
        G = @(x) homotopy.newton(fun,x,var0);
    elseif Opt.homotopy.f2
        % nothing
    else
        error('unknown homotopy-type');
    end
    %
    %% Residual
    %
    if Opt.homotopy.fix || Opt.homotopy.newton
        fH = @(x,l) homotopy.merge(G,fun,x,l,Opt);
        l0 = 0;
        ls = l0;
        le = 1;
        ds0 = 0.01;
    elseif Opt.homotopy.f2
        fH = @(x,l) homotopy.f2(fun,x,l);
        [var0,l0] = homotopy.f2_LM_init(fun,var0);
        if l0>1
            ls = inf;
            le = 1;
            ds0 = abs(le-l0)/10;
        else
            ls = l0;
            le = -inf;
            ds0 = 10^-3;
        end
    else
        error('unknown homotopy-type');
    end
    lt = 1;
    %
    %% Path continuation
    %
    Opt_c = aux.setoff(Opt,'homotopy');
    [var_all,l_all,exitflag] = continuation(fH,var0,ls,le,ds0,'Opt',Opt_c,'l_0',l0,'l_target',lt);
    ind0 = length(l_all)-2;
    [~,ind] = min(l_all(ind0:end)-1);
    ind = ind0-1+ind;
    var_h = var_all(:,ind);
    %
    %% Display
    %
    if exitflag>=0
        aux.print_line(Opt,'success\n');
    else
        aux.print_line(Opt,'failed\n');
    end
    %
end