%% path continuation - draw_corrector
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   25.11.2020 - Alwin Förster
%
function [lco,vco] = draw_corrector( var_all, l_all, dsim1, Opt )
    %% initialize
    nres = 50;
    nv = numel(var_all(:,1));
    lco = cell(1,2);
    vco = cell(1,2);
    lco{1} = NaN(nv,nres);
    vco{1} = NaN(nv,nres);
    lco{2} = NaN(nv,nres);
    vco{2} = NaN(nv,nres);
    %% calc. curves
    if Opt.corrector.sphere
        %% sphere
        %
        if numel(l_all)==1
            lc = l_all;
            vc = var_all;
            rpro = dsim1*ones(nv,1);
        else
            lc = l_all(end-1);
            vc = var_all(:,end-1);
            rpro = sqrt((var_all(:,end)-var_all(:,end-1)).^2+(l_all(end)-l_all(end-1)).^2);
        end
        lco{1} = lc+rpro*cos(linspace(0,2*pi,nres));
        vco{1} = vc+rpro*sin(linspace(0,2*pi,nres));
        dsc = dsim1*ones(nv,1);
        lco{2} = lc+dsc*cos(linspace(0,2*pi,nres));
        vco{2} = vc+dsc*sin(linspace(0,2*pi,nres));
        %
    elseif Opt.corrector.linear
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
            lco{1} = ones(nv,1)*(l_all(end-1)+sec(end)*linspace(0,1,nres));
            vco{1} = var_all(:,end-1)+sec(1:end-1)*linspace(0,1,nres);
            lco{2} = ones(nv,1)*(l_all(end-1)+sec(end)+ort(end)*linspace(0,1,nres));
            vco{2} = var_all(:,end-1)+sec(1:end-1)+ort(1:end-1)*linspace(0,1,nres);
        end
        %
    elseif Opt.corrector.ellipsoid
        %% ellipsoid
        %
%         TODO
        %
    elseif Opt.corrector.ellipsoid2
        %% ellipsoid2
        %
%         TODO
        %
    elseif Opt.corrector.unique
        %% unique
        %
        if numel(l_all)==1
            lc = l_all;
            vc = var_all;
        else
            lc = l_all(end-1);
            vc = var_all(:,end-1);
        end
        dsc = dsim1*sign(Opt.direction(end))*ones(nv,1);
        lco{2} = lc+dsc*ones(1,nres);
        vco{2} = vc+(2*max(abs(var_all).').')*linspace(-1,+1,nres);
        %
    else
        error('corrector-method can not be drawn');
    end
end