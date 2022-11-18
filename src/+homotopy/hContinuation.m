%% path continuation - homotopy.continuation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [ varH, exitflag, varAll, lAll ] = hContinuation( fun, var0, Opt )
    %% Display:
    %
    aux.printLine(Opt,'------> starting homotopy...');
    %
    %% Homotopy-function
    %
    if Opt.homotopy.fix
        G = @(x) homotopy.fix(x,var0);
    elseif Opt.homotopy.newton
        G = @(x) homotopy.newton(fun,x,fun(var0),Opt);
    elseif Opt.homotopy.fixnt
        G = @(x) homotopy.fixNt(fun,x,fun(var0),var0,Opt);
    elseif Opt.homotopy.f2
        % nothing
    elseif Opt.homotopy.squared
        G = @(x) homotopy.squared(x,var0);
    else
        error('unknown homotopy-type');
    end
    %
    %% Residual
    %
    if Opt.homotopy.f2
        fH = @(x,l) homotopy.f2(fun,x,l,Opt);
        [var0,l0] = homotopy.f2LMInit(fun,var0);
        if l0>1
            ls = 2*l0;
        else
            ls = 0;
        end
        le = 1;
        ds0 = abs(le-l0)/100;
    else
        fH = @(x,l) homotopy.merge(G,fun,x,l,Opt);
        l0 = 0;
        ls = -1;
        le = 1;
        ds0 = 0.01;
    end
    lt = 1;
    %
    %% Path continuation
    %
    OptC = aux.setoff(Opt,'homotopy');
    [varAll,lAll,exitflag] = continuation(fH,var0,ls,le,ds0,'Opt',OptC,'l0',l0,'lTarget',lt);
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
        aux.printLine(Opt,'--> homotopy success\n');
    else
        aux.printLine(Opt,'--> homotopy failed\n');
    end
    %
end