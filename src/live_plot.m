%% path continuation - live_plot
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   30.09.2020 - Tido Kubatschek
%
function [pl] = live_plot(Opt, nv, l_all, var_all, pl, bif_flag, bif)
    if length(l_all) == 1
        close all;
        figure;
        clf;
        pl = plot(l_all,var_all,'.-','LineWidth',2);
        grid on;
        xlabel('$\lambda$','interpreter','latex');
        ylabel('$v_{i}$','interpreter','latex');
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
               plot(l_all(bif(1,bif(2,:)==0)),var_all(:,bif(1,bif(2,:)==0)),'bo','LineWidth',2);
               plot(l_all(bif(1,bif(2,:)==1)),var_all(:,bif(1,bif(2,:)==1)),'rx','LineWidth',2);
               hold off;
           end
        end
        drawnow;
    end
end

