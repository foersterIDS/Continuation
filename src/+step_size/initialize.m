%% path continuation - step_size.initialize
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   18.01.2022 - Tido Kubatschek
%
function [Opt,ds0,Stepsize_options] = initialize(Opt,var0,l_start,l_end,ds0)
    %% check chosen stepsize options
    %
    if Opt.check_step_size_options
        %
        if ds0 < Opt.ds_min % ds0 must not be lower than ds_min
            error('ds0 cannot be smaller than ds_min. Consider adapting one of them.');
        end
        %
        if ds0 > Opt.ds_max % ds0 must not be larger than ds_max
            error('ds0 cannot be larger than ds_max. Consider adapting one of them.');
        end
        %
        if Opt.ds_min < 10*Opt.solver_tol % ds_min cannot be lower than tolerance of solver
            Opt.ds_min = 10*Opt.solver_tol;
            aux.print_line(Opt,'--> ds_min has to be at least 10*solver_tol = %.2e. ds_min has been adapted.\n',Opt.ds_min);
        end
        %
        if ds0 > (sqrt(sum(var0.^2) + (l_end-l_start)^2)/10) % ds_max should not be larger than mag of points
            % get order of magnitude
            n_mag = floor(log10(sqrt(sum(var0.^2) + (l_end-l_start)^2)));
            % adapt ds0
            ds0 = 10^(n_mag-1);
            % order of magnitude must be at least 1 lower
            aux.print_line(Opt,'--> ds0 is too large! It must not be greater than 1.00e%i. ds0 has been adapted.\n',n_mag-1);
        end
        %
    end
    %
    %% Tell Stepsize_options which values must be calculated
    
    if Opt.step_size_control.multiplicative || Opt.step_size_control.multiplicative_alt
        %% check weigths
        %
        weights = Opt.weights_multiplicative;
        Stepsize_options = struct('speed_of_continuation',weights(3) ~= 0,...
                                  'predictor', weights(5) ~= 0,...
                                  'rate_of_contraction',weights(4) ~= 0,...
                                  'iterations', true);
    elseif Opt.step_size_control.error || Opt.step_size_control.error_alt
        %% check weigths
        %
        weights = Opt.weights_error;
        Stepsize_options = struct('speed_of_continuation',weights(3) ~= 0,...
                                  'predictor', weights(5) ~= 0,...
                                  'rate_of_contraction',weights(4) ~= 0,...
                                  'iterations', true);
    elseif Opt.step_size_control.contraction
        Stepsize_options = struct('speed_of_continuation',false,...
                                  'predictor',false,...
                                  'rate_of_contraction',true,...
                                  'iterations', true);
    elseif Opt.step_size_control.fix
        Stepsize_options = struct('speed_of_continuation',false,...
                                  'predictor',false,...
                                  'rate_of_contraction',false,...
                                  'iterations', true);
        Opt.ds_max = ds0;
        Opt.ds_min = ds0;
    else
        Stepsize_options = struct('speed_of_continuation',false,...
                                  'predictor',false,...
                                  'rate_of_contraction',false,...
                                  'iterations', true);
    end
    %
end