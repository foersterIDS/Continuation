%% path continuation - homotopy
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [ xr, exitflag ] = homotopy( R, x0, Opt )
    %% Display:
    %
    if Opt.display
        fprintf('------> starting homotopy...');
    end
    %
    %% Homotopy-function
    %
    if Opt.homotopy.fix
        G = @(x) homotopy_fix(x,x0);
    elseif Opt.homotopy.newton
        G = @(x) homotopy_newton(R,x,x0);
    elseif Opt.homotopy.f2
        
    else
        error('unknown homotopy-type');
    end
    %
    %% Residual
    %
    if Opt.homotopy.fix || Opt.homotopy.newton
        fH = @(x,l) homotopy_merge(G,R,x,l,Opt);
        l0 = 0;
        le = 1;
        ds0 = 0.1;
    elseif Opt.homotopy.f2
        fH = @(x,l) homotopy_f2(R,x,l);
        [x0,l0] = homotopy_f2_LM_init(R,x0);
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
    [xs,ll,exitflag] = continuation(fH,x0,ls,le,ds0,'homotopy','off','display','off','l_0',l0,'l_target',lt);
    ind0 = length(ll)-2;
    [~,ind] = min(ll(ind0:end)-1);
    ind = ind0-1+ind;
    xr = xs(:,ind);
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