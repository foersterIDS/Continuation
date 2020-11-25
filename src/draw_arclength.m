%% path continuation - draw_arclength
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   25.11.2020 - Alwin Förster
%
function [lal,val] = draw_arclength( var_all, l_all, ds, dsim1, Opt )
    %% initialize
    nres = 100;
    nv = numel(var_all(:,1));
    lal = cell(1,2);
    val = cell(1,2);
    lal{1} = NaN(nv,nres);
    val{1} = NaN(nv,nres);
    lal{2} = NaN(nv,nres);
    val{2} = NaN(nv,nres);
    %% calc. curves
    if Opt.arclength.sphere
        %% sphere
        %
        if numel(l_all)==1
            lc = l_all;
            vc = var_all;
            dsc = ds*ones(nv,1);
        else
            lc = l_all(end-1);
            vc = var_all(:,end-1);
            dsc = sqrt((var_all(:,end)-var_all(:,end-1)).^2+(l_all(end)-l_all(end-1))^2);
        end
        lal{2} = lc+dsc*cos(linspace(0,2*pi,nres));
        val{2} = vc+dsc*sin(linspace(0,2*pi,nres));
        %
    elseif Opt.arclength.linear
        %% linear
        %
        if numel(l_all)>1
            if numel(l_all)==2
                if numel(Opt.direction)==1
                    sec = [zeros(nv,1);1];
                else
                    sec = Opt.direction;
                end
            else
                sec = [var_all(:,end-1);l_all(end-1)]-[var_all(:,end-2);l_all(end-2)];
            end
            sec = dsim1/norm(sec)*sec;
            ort = ([var_all(:,end);l_all(end)]-[var_all(:,end-1);l_all(end-1)])-sec;
            lal{1} = ones(nv,1)*(l_all(end-1)+sec(end)*linspace(0,1,nres));
            val{1} = var_all(:,end-1)+sec(1:end-1)*linspace(0,1,nres);
            lal{2} = ones(nv,1)*(l_all(end-1)+sec(end)+ort(end)*linspace(0,1,nres));
            val{2} = var_all(:,end-1)+sec(1:end-1)+ort(1:end-1)*linspace(0,1,nres);
        end
        %
    elseif Opt.arclength.ellipsoid
        %% ellipsoid
        %
        error('TODO');
        %
    elseif Opt.arclength.ellipsoid2
        %% ellipsoid2
        %
        error('TODO');
        %
    elseif Opt.arclength.unique
        %% unique
        %
        error('TODO');
        %
    else
        error('arclength-method can not be drawn');
    end
end