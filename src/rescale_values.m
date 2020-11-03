%% path continuation - rescale_values
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   31.10.2020 - Tido Kubatschek
%
%     % Opt needs new entries:
%     % Es gibt jetzt substruct Opt_scaling --> Opt.scaling mit dem Eintrag
%     % 'dynamicDscale',false --> auf true setzen wenn das hier
%     % funktioniert

%     % Dscale, Dsacle0 on of which cannot be filled directly by the user
%     % --> Kann ruhig beides bei Opt rein. Ich mache mir mal gedanken wie man
%     % da Einträge sperrt.
%     
function [x_scaled, Opt, jacobian_scaled] = rescale_values(Opt,var_all,l_all,jacobian)
    %% RESCALE (dynamically adjust scaling, if requested)
    if ison(Opt.scaling)
        x_all = [var_all; l_all];
        if length(l_all) == 1
            Opt.Dscale0 = Opt.Dscale;
        end
        %% Adjust variable scaling values (diagonal elements of diagonal
        % scaling matrix) so that scaled variables have value ~1. But avoid
        % weighing very small values too much by setting as minimum dynamic
        % scaling the initial one (Dscale0).
        
        Dscaleold = Opt.Dscale;
        Opt.Dscale(1:end-1) = max(abs(x_all(1:end-1,end).*Dscaleold(1:end-1)),...
            Opt.Dscale0(1:end-1));
        
        %% Update scaled variable values, Jacobian, reference tangent
        x_scaled      = x_all(:,end).*(Dscaleold./Opt.Dscale);
%         Xref    = Xref.*(Dscaleold./Otp.Dscale);
%         Xold    = Xold.*(Dscaleold./Opt.Dscale);
%         jacobian_scaled = jacobian*diag(1./Dscaleold);
%         zref    = zref.*Dscaleold;
        
    end
end

