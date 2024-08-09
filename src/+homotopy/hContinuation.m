%% path continuation - homotopy.continuation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [ varH, exitflag, varAll, lAll ] = hContinuation( fun, var0, oih )
    %% Display:
    %
    aux.printLine(oih,'------> starting homotopy...');
    %
    %% Homotopy-function
    %
    if oih.opt.homotopy.fix
        G = @(x) homotopy.fix(x,var0);
    elseif oih.opt.homotopy.newton
        G = @(x) homotopy.newton(fun,x,fun(var0),oih);
    elseif oih.opt.homotopy.fixnt
        G = @(x) homotopy.fixNt(fun,x,fun(var0),var0,oih);
    elseif oih.opt.homotopy.f2
        % nothing
    elseif oih.opt.homotopy.squared
        G = @(x) homotopy.squared(x,var0);
    else
        error('unknown homotopy-type');
    end
    %
    %% Residual
    %
    if oih.opt.homotopy.f2
        fH = @(x,l) homotopy.f2(fun,x,l,oih);
        [var0,l0] = homotopy.f2LMInit(fun,var0);
        if l0>1
            ls = 2*l0;
        else
            ls = 0;
        end
        le = 1;
        ds0 = abs(le-l0)/100;
    else
        fH = @(x,l) homotopy.merge(G,fun,x,l,oih);
        l0 = 0;
        ls = -1;
        le = 1;
        ds0 = 0.01;
    end
    lt = 1;
    %
    %% Path continuation
    %
    oih.opt = aux.setoff(oih.opt,'homotopy');
    [varAll,lAll,exitflag] = continuation(fH,var0,ls,le,ds0,'Opt',OptC,'l0',l0,'lTarget',lt);
    oih.opt = aux.seton(oih.opt,'homotopy');
    ind0 = numel(lAll)-2;
    if ind0>=1
        [~,ind] = min(lAll(ind0:end)-1);
        ind = ind0-1+ind;
    elseif numel(lAll)==0
        ind = [];
    else
        ind  = numel(lAll);
    end
    varH = varAll(:,ind);
    %
    %% Display
    %
    if exitflag>0
        aux.printLine(oih,'--> homotopy success\n');
    else
        aux.printLine(oih,'--> homotopy failed\n');
    end
    %
end