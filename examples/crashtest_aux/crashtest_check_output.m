%% path continuation - validate_result
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   23.10.2020 - Alwin Förster
%
function [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,info,probinfo,probcounter)
    [~,ind] = find(isnan(ls));
    if (exitflag<=0) || (abs(ls(end)-lame)/abs(lame-lams)>10^-3 && abs(ls(end)-lams)/abs(lams-lame)>10^-3)
        probinfo = [probinfo,info];
        probcounter = probcounter+1;
        if abs(ls(end)-lame)/abs(lame-lams)>10^-3 && abs(ls(end)-lams)/abs(lams-lame)>10^-3
            probinfo = [probinfo,sprintf('--> continuation has not reached the target: %.2f%% complete\n',(ls(end)-lams)/(lame-lams)*100)];
        end
        if exitflag<=0
            probinfo = [probinfo,sprintf('--> continuation exitflag: %d\n',exitflag)];
        end
    end
end