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
function [var_all,l_all,exitflag,bif] = continuation(fun,var0,l_start,l_end,ds0,varargin)
    %% initialize
    %
    exitflag = -1;
    Opt = continuation_input(varargin,fun,var0,l_start,l_end);
    solver = continuation_solver(Opt);
    res_arle = residual_arclength(Opt);
    ds = ds0;
    do_deflate = false;
    do_homotopy = false;
    error_counter = 0;
    loop_counter = 0;
    if Opt.display
        fprintf('Starting path continuation...\n');
        t_display = tic;
    end
    %
    %% find initial solution
    %
    residual_initial = @(v) residual_fixed_value(fun,v,Opt.l_0,Opt);
    [var_all,~,initial_exitflag] = solver(residual_initial,var0);
    if initial_exitflag>0
        l_all = Opt.l_0;
        do_continuation = true;
        if Opt.display
            fprintf('Initial solution at l = %.2e\n',Opt.l_0);
        end
    else
        var_all = [];
        l_all = [];
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
        residual = @(x) merge_residuals(fun,res_arle,x,[var_all;l_all],ds,Opt);
        if do_deflate
            try
                residual = @(x) deflation(residual,x_deflation,x,Opt);
            catch
                error('Error occured during deflation.');
            end
        end
        [v_predictor,l_predictor] = predictor(var_all,l_all,ds,Opt);
        x_predictor = [v_predictor;l_predictor];
        %
        %% solve
        %
        if do_homotopy
            try
                x_last_step = [var_all(:,end);l_all(end)];
                [x_solution,solver_exitflag] = homotopy(residual,x_last_step,Opt);
            catch
                x_solution = NaN(size(x_predictor));
                solver_exitflag = -2;
            end
        else
            try
                [x_solution,~,solver_exitflag] = solver(residual,x_predictor);
            catch
                x_solution = NaN(size(x_predictor));
                solver_exitflag = -2;
            end
        end
        %
        %% check result
        %
        val = validate_result(x_solution,var_all,l_all,solver_exitflag,Opt);
        if val
            %% valid result
            var_all = [var_all,x_solution(1:end-1)];
            l_all = [l_all,x_solution(end)];
            do_deflate = false;
            do_homotopy = false;
            error_counter = 0;
        else
            %% invalid result
            error_counter = error_counter+1;
            if Opt.deflation && ~isnan(sum(x_solution(:,end)))
                do_deflate = true;
                x_deflation = x_solution;
            else
                do_deflate = false;
                warning('Hier muss etwas gemacht werden!');
            end
            if ison(Opt.homotopy) && error_counter>=Opt.homotopy_error_counter
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
                fprintf('-----> continued at l = %.2e\t|\tnew arc-length: ds = %.2e\t|\tloop counter = %d\n',l_all(end),ds,loop_counter);
            else
                fprintf('-----> invalid point\t\t\t\t\t\t\t\t\t\t\t\t\t|\tloop counter = %d\n',loop_counter);
            end
        end
        % exit with success:
        if l_all(end)>l_end
            do_continuation = false;
            exitflag = 1;
        end
        % exit with l<l_start:
        if l_all(end)<l_start
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
    %% final disp
    %
    if Opt.display
        fprintf('Time Elapsed: %.3f s\n',toc(t_display));
    end
    %
end