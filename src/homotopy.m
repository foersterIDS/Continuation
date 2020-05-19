%% path continuation - homotopy
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin F�rster
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
        G = @(x) x-x0;
    elseif Opt.homotopy.newton
        G = @(x) R(x)-R(x0);
    else
        error('unknown homotopy-type');
    end
    %
    %% Residual
    %
    fH = @(x,l) (1-l)*G(x)+l*R(x);
    %
    %% Path continuation
    %
    ds0 = 0.1; % TODO
    [xs,~,exitflag] = continuation(fH,x0,0,1,ds0,'homotopy','off','display','off');
    xr = xs(:,end);
    %
    %% Display
    %
    if Opt.display
        if exitflag>0
            fprintf('success\n');
        else
            fprintf('failed\n');
        end
    end
    %
end