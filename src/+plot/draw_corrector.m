%% path continuation - plot.draw_corrector
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   25.11.2020 - Alwin Förster
%
function [lco,vco] = draw_corrector( Path, dsim1, Opt )
    %% initialize
    nres = 50;
    nv = numel(Path.var_all(:,1));
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
        if numel(Path.l_all)==1
            lc = Path.l_all;
            vc = Path.var_all;
            rpro = dsim1*ones(nv,1);
        else
            lc = Path.l_all(end-1);
            vc = Path.var_all(:,end-1);
            rpro = sqrt((Path.var_all(:,end)-Path.var_all(:,end-1)).^2+(Path.l_all(end)-Path.l_all(end-1)).^2);
        end
        lco{1} = lc+rpro*cos(linspace(0,2*pi,nres));
        vco{1} = vc+rpro*sin(linspace(0,2*pi,nres));
        dsc = dsim1*ones(nv,1);
        lco{2} = lc+dsc*cos(linspace(0,2*pi,nres));
        vco{2} = vc+dsc*sin(linspace(0,2*pi,nres));
        %
    elseif Opt.corrector.orthogonal
        %% orthogonal
        %
        if numel(Path.l_all)>1
            if numel(Path.l_all)==2
                if numel(Opt.direction)==1
                    sec = [zeros(nv,1);1];
                else
                    sec = Opt.direction;
                end
            else
                sec = [Path.var_all(:,end-1);Path.l_all(end-1)]-[Path.var_all(:,end-2);Path.l_all(end-2)];
            end
            sec = dsim1/norm(sec)*sec;
            ort = ([Path.var_all(:,end);Path.l_all(end)]-[Path.var_all(:,end-1);Path.l_all(end-1)])-sec;
            lco{1} = ones(nv,1)*(Path.l_all(end-1)+sec(end)*linspace(0,1,nres));
            vco{1} = Path.var_all(:,end-1)+sec(1:end-1)*linspace(0,1,nres);
            lco{2} = ones(nv,1)*(Path.l_all(end-1)+sec(end)+ort(end)*linspace(0,1,nres));
            vco{2} = Path.var_all(:,end-1)+sec(1:end-1)+ort(1:end-1)*linspace(0,1,nres);
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
        if numel(Path.l_all)==1
            lc = Path.l_all;
            vc = Path.var_all;
        else
            lc = Path.l_all(end-1);
            vc = Path.var_all(:,end-1);
        end
        dsc = dsim1*sign(Opt.direction(end))*ones(nv,1);
        lco{2} = lc+dsc*ones(1,nres);
        vco{2} = vc+(2*max(abs(Path.var_all).').')*linspace(-1,+1,nres);
        %
    else
        error('corrector-method can not be drawn');
    end
end