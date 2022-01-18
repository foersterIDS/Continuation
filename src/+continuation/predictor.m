%% path continuation - continuation.predictor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [vp,lp,fun_predictor,sp,ds] = predictor(Path,ds,solver_jacobian,fun,res_corr,Solver,Opt)
    %% get fun_predictor:
    if Opt.predictor.polynomial
        if length(Path.l_all)==1
            fun_predictor = @(s) predictor.initial(Path,s,Opt);
        else
            [nt,nf] = predictor.adaptive(Path,Opt);
            [fpt,Jpt] = predictor.taylor(Path,nt,nf);
            fun_predictor = @(s) aux.fncHndToVal(s,fpt,Jpt);
        end
    elseif Opt.predictor.tangential
        if length(Path.l_all)==1
            fun_predictor = @(s) predictor.initial(Path,s,Opt);
        else
            fun_predictor = @(s) predictor.ode(Path,s,solver_jacobian,fun);
        end
    else
        error('predictor not set or of unknown type');
    end
    %% predictor_solver:
    if Opt.predictor_solver
        xi = [Path.var_all(:,end);Path.l_all(end)];
        if Opt.enforce_ds_max
            %% enforce_ds_max:
            p = 0.95;
            fun_pred = @(s,ds) aux.merge_arle_pred(fun_predictor,res_corr,s,xi,ds);
            fun_ds_max = @(s,ds) heaviside(-min(Opt.ds_max*p-(fun_predictor(s)-xi)))*min(Opt.ds_max*p-(fun_predictor(s)-xi));
            fun_solve = @(spds) [fun_pred(spds(1),spds(2));fun_ds_max(spds(1),spds(2))];
            [spds,~,exitflag] = Solver.num_jac(fun_solve,[ds;ds]);
            if exitflag<=0
                sp = ds;
            else
                sp = spds(1);
                ds = spds(2);
            end
        else
            %% solve arc-length
            fun_solve = @(s) aux.merge_arle_pred(fun_predictor,res_corr,s,xi,ds);
            [sp,~,exitflag] = Solver.predictor(fun_solve,ds);
            if exitflag<=0
                sp = ds;
            end
        end
    else
        sp = ds;
    end
    %% get predictor:
    xip1 = fun_predictor(sp);
    %% correct_predictor:
    if Opt.correct_predictor && numel(Path.l_all)>1
        dxi = [Path.var_all(:,end);Path.l_all(end)]-[Path.var_all(:,end-1);Path.l_all(end-1)];
        dxip1 = xip1-[Path.var_all(:,end);Path.l_all(end)];
        if dot(dxip1,dxi)<0
            dxip1 = -dxip1;
            xip1 = [Path.var_all(:,end);Path.l_all(end)]+dxip1;
        end
    end
    %% make output:
    vp = xip1(1:end-1);
    lp = xip1(end);
end