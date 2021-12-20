%% path continuation - initialize_structs
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   16.09.2021 - Alwin Förster
%
function [Bifurcation,Counter,Do,Info,Path,Plot,Solver] = initialize_structs(var0,l_start,l_end,Opt)
    %% Bifurcation
    %
    Bifurcation = struct('bif',[],...
                         'dirs',{cell(1,2)},...
                         'flag',0);
    %
    %% Counter
    %
    Counter = struct('catch',0,...
                     'catch_old',0,...
                     'error',0,...
                     'loop',0,...
                     'step',0,...
                     'valid_stepback',0);
    %
    %% Do
    %
    Do = struct('continuation',false,...
                'convergeToTarget',false,...
                'deflate',false,...
                'homotopy',false,...
                'loop',false,...
                'remove',false,...
                'stepback',false,...
                'suspend',false);
	%
    %% Info
    %
    Info = struct('exitflag',-1,...
                  'nv',numel(var0),...
                  'l_start',l_start,...
                  'l_end',l_end);
	%
    %% Path
    %
    Path = struct('var_all',[],...
                  'l_all',[],...
                  's_all',[]);
	%
    %% Plot
    %
    if (ison(Opt.plot) && ~Opt.plot.detail) || ~ison(Opt.plot)
        Plot = struct('fig',[],...
                      'pl',[],...
                      'pl_curr',[]);
    else
        Plot = struct('fig',[],...
                      'pl',[],...
                      'pl_it',[],...
                      'pl_det',[],...
                      'pl_s',[],...
                      'pl_cor_assist',[],...
                      'pl_cor',[],...
                      'pl_pre',[],...
                      'pl_curr',[]);
    end
	%
    %% Solver
    %
    [solver,predictor_solver,num_jac_solver,default_solver_output] = continuation_solver(Opt);
    Solver = struct('main',solver,...
                    'predictor',predictor_solver,...
                    'num_jac',num_jac_solver,...
                    'default_output',default_solver_output);
    %
end