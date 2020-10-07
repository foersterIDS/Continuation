%% path continuation - live_plot
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   30.09.2020 - Tido Kubatschek
%
function [pl] = live_plot(Opt, nv, l_start, l_end, l_all, var_all, pl, bif_flag, bif)
    l_lu = [min([l_start,l_end]),max([l_start,l_end])];
    dl = 0.1;
    if length(l_all) == 1
        figure('units', 'normalized', 'position', [0.2,0.3,0.6,0.5]);
        clf;
        pl = plot(l_all,var_all,'.-','LineWidth',2);
        grid on;
        xlabel('$\lambda$','interpreter','latex');
        ylabel('$v_{i}$','interpreter','latex');
        xlim([max([l_lu(1),min(l_all-10^-15)*(1-dl)]),min([l_lu(2),max(l_all+10^-15)*(1+dl)])]);
        drawnow;
    else
        newXData = cell(nv, 1);
        newYData = cell(nv, 1);
        for k = 1:nv
            newXData{k,1} = l_all;
            newYData{k,1} = var_all(k,:);
        end
        set(pl, {'XData'}, newXData, {'YData'},  newYData);
        if ~Opt.unique && ison(Opt.bifurcation) && bif_flag 
           if ~isempty(bif)
               hold on;
               if bif(2,end) == 0
                   plot(l_all(end),var_all(:,end),'bo','LineWidth',2);
               elseif bif(2,end) == 1
                   plot(l_all(end),var_all(:,end),'rx','LineWidth',2);
               end
               hold off;
           end
        end
        xlim([max([l_lu(1),min(l_all-10^-15)*(1-dl)]),min([l_lu(2),max(l_all+10^-15)*(1+dl)])]);
        drawnow;
    end
end

