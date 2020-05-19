%% path continuation - continuation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
%   fun = fun(var,l) != 0
%   lams <= l <= lame
%   ds0: initial stepsize
%
function [vars,ls,exitflag,bif] = continuation(fun,var0,l_start,l_end,ds0,varargin)
    %% initialize
    %
    exitflag = -1;
    Opt = continuation_input(varargin,fun);
    solver = continuation_solver(Opt);
    res_arle = residual_arclength(Opt);
    ds = ds0;
    do_deflate = false;
    do_homotopy = false;
    error_counter = 0;
    loop_counter = 0;
    %
    %% find initial solution
    %
    R = @(v) fun(v,l_start);
    [vars,~,initial_exitflag] = solver(R,var0);
    if initial_exitflag>0
        ls = l_start;
        do_continuation = true;
        if Opt.display
            fprintf('Initial solution at l = %.2e\n',l_start);
        end
    else
        vars = [];
        ls = [];
        bif = [];
        exitflag = -2; % no initial solution found
        do_continuation = false;
        if Opt.display
            fprintf('No initial solution found.\n');
        end
    end
    %
    %% continuation
    %
    while do_continuation
        %% initialize loop
        %
        loop_counter = loop_counter+1;
        %
        %% residual and predictor
        %
        R = @(x) [fun(x(1:end-1),x(end));res_arle(x,[vars;ls],ds)];
        if do_deflate
            try
                R = @(x) deflation(R,xdef,x);
            catch
                error('Error occured during deflation.');
            end
        end
        [vp,lp] = predictor(vars,ls,ds,Opt);
        xp = [vp;lp];
        %
        %% solve
        %
        if do_homotopy
            try
                xi = [vars(:,end);ls(end)];
                [x_sol,solver_exitflag] = homotopy(R,xi,Opt);
            catch
                x_sol = NaN(size(xp));
                solver_exitflag = -2;
            end
        else
            try
                [x_sol,~,solver_exitflag] = solver(R,xp);
            catch
                x_sol = NaN(size(xp));
                solver_exitflag = -2;
            end
        end
        %
        %% check result
        %
        val = validate_result(x_sol,vars,ls,solver_exitflag);
        if val
            %% valid result
            vars = [vars,x_sol(1:end-1)];
            ls = [ls,x_sol(end)];
            do_deflate = false;
            do_homotopy = false;
            error_counter = 0;
        else
            %% invalid result
            error_counter = error_counter+1;
            if Opt.deflation && ~isnan(sum(x_sol(:,end)))
                do_deflate = true;
                xdef = x_sol;
            else
                do_deflate = false;
                warning('Hier muss etwas gemacht werden!');
            end
            if sum(structfun(@(x) x,Opt.homotopy)) && error_counter>=Opt.homotopy_error_counter
                do_homotopy = true;
            else
                do_homotopy = false;
            end
        end
        %
        %% Bifurcations
        %
        bif = [];
        % TODO: Bifurkationen erkennen und exakten Punkt ermitteln
        %
        %% adjust arc-length
        %
        ds = arclength(ds,ds0,error_counter);
        %
        %% end loop
        %
        if Opt.display
            if val
                fprintf('-----> continued at l = %.2e\t|\tnew arc-length: ds = %.2e\t|\tloop counter = %d\n',ls(end),ds,loop_counter);
            else
                fprintf('-----> invalid point\t\t\t\t\t\t\t\t\t\t\t\t\t|\tloop counter = %d\n',loop_counter);
            end
        end
        % exit with success:
        if ls(end)>l_end
            do_continuation = false;
            exitflag = 1;
        end
        % exit with l<l_start:
        if ls(end)<l_start
            do_continuation = false;
            exitflag = 0;
        end
        % exit without complete results:
        if error_counter>=Opt.max_error_counter
            exitflag = -1;
            warning('No valid result could be found for the last %d attempts.',Opt.max_error_counter);
            do_continuation = false;
        end
        %
    end
    %
end