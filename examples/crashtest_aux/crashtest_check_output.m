%% path continuation - validate_result
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   23.10.2020 - Alwin Förster
%
function [probinfo] = crashtest_check_output(vs,ls,exitflag,lams,lame,info,probinfo)
    if (exitflag<=0) || (abs(ls(end)-lame)>10^-8)
        probinfo = [probinfo,info];
        if abs(ls(end)-lame)>10^-8
            probinfo = [probinfo,sprintf('--> continuation has not reached the target: %.2f%% complete\n',(ls(end)-lams)/(lame-lams)*100)];
        end
        if exitflag<=0
            probinfo = [probinfo,sprintf('--> continuation exitflag: %d\n',exitflag)];
        end
    end
end