%% path continuation - homotopy_continuation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [ var_h, exitflag ] = homotopy_continuation( fun, var0, Opt )
    %% Display:
    %
    if Opt.display
        fprintf('------> starting homotopy...');
    end
    %
    %% Homotopy-function
    %
    if Opt.homotopy.fix
        G = @(x) homotopy_fix(x,var0);
    elseif Opt.homotopy.newton
        G = @(x) homotopy_newton(fun,x,var0);
    elseif Opt.homotopy.f2
        % nothing
    else
        error('unknown homotopy-type');
    end
    %
    %% Residual
    %
    if Opt.homotopy.fix || Opt.homotopy.newton
        fH = @(x,l) homotopy_merge(G,fun,x,l,Opt);
        l0 = 0;
        ls = l0;
        le = 1;
        ds0 = 0.1;
    elseif Opt.homotopy.f2
        fH = @(x,l) homotopy_f2(fun,x,l);
        [var0,l0] = homotopy_f2_LM_init(fun,var0);
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
    Opt_c = Opt;
    [var_all,l_all,exitflag] = continuation(fH,var0,ls,le,ds0,'Opt',Opt_c,'homotopy','off','l_0',l0,'l_target',lt);
    ind0 = length(l_all)-2;
    [~,ind] = min(l_all(ind0:end)-1);
    ind = ind0-1+ind;
    var_h = var_all(:,ind);
    %
    %% Display
    %
    if Opt.display
        if exitflag>=0
            fprintf('success\n');
        else
            fprintf('failed\n');
        end
    end
    %
end