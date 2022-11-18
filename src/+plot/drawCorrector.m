%% path continuation - plot.drawCorrector
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   25.11.2020 - Alwin Förster
%
function [lco,vco] = drawCorrector( Path, dsim1, Opt )
    %% initialize
    nres = 50;
    nv = numel(Path.varAll(:,1));
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
        if numel(Path.lAll)==1
            lc = Path.lAll;
            vc = Path.varAll;
            rpro = dsim1*ones(nv,1);
        else
            lc = Path.lAll(end-1);
            vc = Path.varAll(:,end-1);
            rpro = sqrt((Path.varAll(:,end)-Path.varAll(:,end-1)).^2+(Path.lAll(end)-Path.lAll(end-1)).^2);
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
        if numel(Path.lAll)>1
            if numel(Path.lAll)==2
                if numel(Opt.direction)==1
                    sec = [zeros(nv,1);1];
                else
                    sec = Opt.direction;
                end
            else
                sec = [Path.varAll(:,end-1);Path.lAll(end-1)]-[Path.varAll(:,end-2);Path.lAll(end-2)];
            end
            sec = dsim1/norm(sec)*sec;
            ort = ([Path.varAll(:,end);Path.lAll(end)]-[Path.varAll(:,end-1);Path.lAll(end-1)])-sec;
            lco{1} = ones(nv,1)*(Path.lAll(end-1)+sec(end)*linspace(0,1,nres));
            vco{1} = Path.varAll(:,end-1)+sec(1:end-1)*linspace(0,1,nres);
            lco{2} = ones(nv,1)*(Path.lAll(end-1)+sec(end)+ort(end)*linspace(0,1,nres));
            vco{2} = Path.varAll(:,end-1)+sec(1:end-1)+ort(1:end-1)*linspace(0,1,nres);
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
        if numel(Path.lAll)==1
            lc = Path.lAll;
            vc = Path.varAll;
        else
            lc = Path.lAll(end-1);
            vc = Path.varAll(:,end-1);
        end
        dsc = dsim1*sign(Opt.direction(end))*ones(nv,1);
        lco{2} = lc+dsc*ones(1,nres);
        vco{2} = vc+(2*max(abs(Path.varAll).').')*linspace(-1,+1,nres);
        %
    else
        error('corrector-method can not be drawn');
    end
end