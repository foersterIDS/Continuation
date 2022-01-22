%% path continuation - fun_solver.basic_newton
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin F�rster
%   16.09.2020 - Tido Kubatschek 
%
function [x,fval,exitflag,output,jacobian] = basic_newton(fun,x0,dscale,output_flag,Opt)
    max_step = Opt.solver_max_iterations;
    if output_flag
        global solver_stepsizes;
        solver_stepsizes = [0, 0];
    end
    n_steps = 0;
    exitflag = -1;
    xi = x0;
    xim1 = x0;
    xi_sc = x0./dscale;
    for i=1:max_step
        if Opt.jacobian
            [fi,Ji] = fun(xi);
        else
            fi = fun(xi);
            Ji = aux.numeric_jacobian(fun,xi,'central_value',fi,'diffquot',Opt.diffquot);
        end
        Ji_sc = Ji*diag(dscale);
        nf = length(fi);
        xip1_sc = xi_sc-Ji_sc(1:nf,1:nf)\fi;
        %% end loop
        xi_sc = xip1_sc;
        xi = xip1_sc.*dscale;
        n_steps = n_steps+1;
        abs_fi = sqrt(fi'*fi);
        abs_dxi = norm(xi-xim1)/norm(xi);
        if output_flag
            solver_stepsizes = [solver_stepsizes; i, norm(xi-xim1)];
        end
        if abs_fi<Opt.solver_tol
            exitflag = 2;
            break;
        elseif abs_dxi<Opt.solver_tol
            exitflag = 1;
            break;
        end
        xim1 = xi;
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
                    Ji = aux.numeric_jacobian(fun,xi,'central_value',fi,'diffquot',Opt.diffquot);
                end
                jacobian = Ji;
            end
        end
    end
end