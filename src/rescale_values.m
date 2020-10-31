%% path continuation - rescale_values
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   31.10.2020 - Tido Kubatschek
%
function [] = rescale_values(Opt)
%     % Opt needs new entries:
%     % dynamicDscale
%     % Dscale, Dsacle0 on of which cannot be filled directly by the user
%     
%     %% RESCALE (dynamically adjust scaling, if requested)
%     if Opt.dynamicDscale
%        
%         %% Adjust variable scaling values (diagonal elements of diagonal
%         % scaling matrix) so that scaled variables have value ~1. But avoid
%         % weighing very small values too much by setting as minimum dynamic
%         % scaling the initial one (Dscale0).
%         Dscaleold = Opt.Dscale;
%         Opt.Dscale(1:end-1) = max(abs(X0(1:end-1).*Dscaleold(1:end-1)),...
%             Opt.Dscale0(1:end-1));
%         
%         %% Update scaled variable values, Jacobian, reference tangent
%         Xref    = Xref.*(Dscaleold./Otp.Dscale);
%         X0      = X0.*(Dscaleold./Opt.Dscale);
%         Xold    = Xold.*(Dscaleold./Opt.Dscale);
%         if Sopt.flag
%             J       = J*diag(1./Dscaleold);
%             zref    = zref.*Dscaleold;
%         end
%     end
end

