%% path continuation - step_size.contraction
%  Adjusts stepsize by the rate of contraction of the solver which is
%  compared to an optimale rate.
%
%
%   Inputs:
%       Solver.output -- contains information of solver, such as the 
%                        rate of contraction.
%       Opt           -- contains user inputs, such as the optimal contraction
%                        rate specified in 'optimal_contraction_rate'.
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('step_size.control')">other stepsize adaption methods</a>.
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.01.2022 - Tido Kubatschek
%
function [xi] = contraction(Solver,Opt)
    if ~isempty(Solver.output.rate_of_contraction)
        % get rates of contraction
        %
        % current rate
        %
        current_rate = Solver.output.rate_of_contraction(end);
        %
        % if there is more than one rate accessible also get previous rate
        % and calculate relative difference
        %
        if length(Solver.output.rate_of_contraction) > 1
            previous_rate = Solver.output.rate_of_contraction(end-1);
            rel_difference = abs(current_rate - previous_rate) / abs(previous_rate);
        else
            rel_difference = inf;
        end
        % 
        % optimal rate
        %
        q_opt = Opt.optimal_contraction_rate;
        %
        if current_rate > q_opt && rel_difference > 0.01
            % only decrease if there was a difference
            xi = 0.5;
        elseif current_rate <= q_opt/4 || (current_rate <= q_opt/2 && rel_difference < 0.01)
            % increase if rate is small or there is just a small difference
            xi = 1.2;
        else
            xi = 1;
        end
    else
        xi = 1;
    end
end