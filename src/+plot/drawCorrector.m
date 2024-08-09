%% path continuation - plot.drawCorrector
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   25.11.2020 - Alwin Förster
%
function [lco,vco] = drawCorrector( oih, dsim1 )
    %% initialize
    nres = 50;
    nv = numel(oih.path.varAll(:,1));
    lco = cell(1,2);
    vco = cell(1,2);
    lco{1} = NaN(nv,nres);
    vco{1} = NaN(nv,nres);
    lco{2} = NaN(nv,nres);
    vco{2} = NaN(nv,nres);
    %% calc. curves
    if oih.opt.corrector.sphere
        %% sphere
        %
        if oih.path.nAll==1
            lc = oih.path.lAll;
            vc = oih.path.varAll;
            rpro = dsim1*ones(nv,1);
        else
            lc = oih.path.lAll(end-1);
            vc = oih.path.varAll(:,end-1);
            rpro = sqrt((oih.path.varAll(:,end)-oih.path.varAll(:,end-1)).^2+(oih.path.lAll(end)-oih.path.lAll(end-1)).^2);
        end
        lco{1} = lc+rpro*cos(linspace(0,2*pi,nres));
        vco{1} = vc+rpro*sin(linspace(0,2*pi,nres));
        dsc = dsim1*ones(nv,1);
        lco{2} = lc+dsc*cos(linspace(0,2*pi,nres));
        vco{2} = vc+dsc*sin(linspace(0,2*pi,nres));
        %
    elseif oih.opt.corrector.orthogonal
        %% orthogonal
        %
        if oih.path.nAll>1
            if oih.path.nAll==2
                if numel(oih.opt.direction)==1
                    sec = [zeros(nv,1);1];
                else
                    sec = oih.opt.direction;
                end
            else
                sec = [oih.path.varAll(:,end-1);oih.path.lAll(end-1)]-[oih.path.varAll(:,end-2);oih.path.lAll(end-2)];
            end
            sec = dsim1/norm(sec)*sec;
            ort = ([oih.path.varAll(:,end);oih.path.lAll(end)]-[oih.path.varAll(:,end-1);oih.path.lAll(end-1)])-sec;
            lco{1} = ones(nv,1)*(oih.path.lAll(end-1)+sec(end)*linspace(0,1,nres));
            vco{1} = oih.path.varAll(:,end-1)+sec(1:end-1)*linspace(0,1,nres);
            lco{2} = ones(nv,1)*(oih.path.lAll(end-1)+sec(end)+ort(end)*linspace(0,1,nres));
            vco{2} = oih.path.varAll(:,end-1)+sec(1:end-1)+ort(1:end-1)*linspace(0,1,nres);
        end
        %
    elseif oih.opt.corrector.ellipsoid
        %% ellipsoid
        %
%         TODO
        %
    elseif oih.opt.corrector.ellipsoid2
        %% ellipsoid2
        %
%         TODO
        %
    elseif oih.opt.corrector.unique
        %% unique
        %
        if oih.path.nAll==1
            lc = oih.path.lAll;
            vc = oih.path.varAll;
        else
            lc = oih.path.lAll(end-1);
            vc = oih.path.varAll(:,end-1);
        end
        dsc = dsim1*sign(oih.opt.direction(end))*ones(nv,1);
        lco{2} = lc+dsc*ones(1,nres);
        vco{2} = vc+(2*max(abs(oih.path.varAll).').')*linspace(-1,+1,nres);
        %
    else
        error('corrector-method can not be drawn');
    end
end