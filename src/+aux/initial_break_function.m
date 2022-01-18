%% path continuation - aux.initial_break_function
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   21.02.2021 - Alwin Förster
%
function [do_break,break_fun_out] = initial_break_function(fun,jac,v,l,break_fun_out)
    do_break = false;
    break_fun_out = [break_fun_out,numel(break_fun_out)+1];
end