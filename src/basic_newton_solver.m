%% path continuation - basic_newton_solver
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%   16.09.2020 - Tido Kubatschek 
%
function [x,fval,exitflag,output,jacobian] = basic_newton_solver(fun,x0,Opt)
    tol = 10^-8;
    max_step = 50;
    n_steps = 0;
    exitflag = -1;
    xi = x0;
    for i=1:max_step
        if Opt.jacobian
            [fi,Ji] = fun(xi);
        else
            fi = fun(xi);
            Ji = numeric_jacobian(fun,xi,'fx',fi);
        end
        nf = length(fi);
        xip1 = xi-Ji(1:nf,1:nf)\fi;
        abs_fi = sqrt(fi'*fi);
        if abs_fi<tol
            exitflag = 1;
            break;
        end
        %% end loop
        xi = xip1;
        n_steps = n_steps+1;
    end
    %% make output
    x = xi;
    if n_steps==max_step
        % TODO: hier muss es noch mehr geben
        exitflag = 0;
    end
    if nargout>1
        if Opt.jacobian
            [fi,Ji] = fun(xi);
        else
            fi = fun(xi);
        end
        fval = fi;
        if nargout>3
            output = struct('iterations',n_steps,...
                'algorithm','basic-newton',...
                'tolerance',abs_fi);
            if nargout>4
                if ~Opt.jacobian
                    Ji = numeric_jacobian(fun,xi,'fx',fi);
                end
                jacobian = Ji;
            end
        end
    end
end