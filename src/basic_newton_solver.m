%% path continuation - basic_newton_solver
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%
function [x,fval,exitflag,output,jacobian] = basic_newton_solver(fun,x0)
    is_jacobian = abs(nargout(fun))==2;
    tol = 10^-8;
    max_step = 50;
    n_steps = 0;
    exitflag = -1;
    xi = x0;
    for i=1:max_step
        if is_jacobian
            [fi,Ji] = fun(xi);
        else
            fi = fun(xi);
            Ji = numeric_jacobian(fun,xi,fi);
        end
        xip1 = xi-Ji\fi;
        toli = sqrt((xip1-xi)'*(xip1-xi));
        if toli<tol
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
        if is_jacobian
            [fi,Ji] = fun(xi);
        else
            fi = fun(xi);
        end
        fval = fi;
    	if nargout>3
            output = struct('iterations',n_steps,...
                            'algorithm','basic-newton',...
                            'tolerance',toli);
            if nargout>4
                if ~is_jacobian
                    Ji = numeric_jacobian(fun,xi,fi);
                end
                jacobian = Ji;
            end
        end
    end
end