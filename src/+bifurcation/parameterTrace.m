%% path continuation - bifurcation.parameterTrace
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   31.03.2022 - Alwin Förster
%
function [Path] = parameterTrace(Opt,Path,Bifurcation,Info,fun)
    %% settings for dpa:
    bifTrace = Bifurcation.bif(1,:);
    nbif = numel(bifTrace);
    OptTrace = Opt;
    OptTrace = aux.setoff(OptTrace,'bifurcation');
    OptTrace.direction = [zeros(size(OptTrace.direction));sign(Opt.gTarget-Opt.g0)];
    OptTrace.l0 = Opt.g0;
    OptTrace.lTarget = Opt.gTarget;
    OptTrace.dpaGammaVar = true;
    PathBifs = struct('varAll',[],'lAll',[],'sAll',[]);
    %% start dpa:
    for ii=1:nbif
        i = bifTrace(ii);
        sc = Bifurcation.scaling(ii);
        resDpa = @(v,l,g) dpa.resBif(fun,[v;l],g,Opt,sc);
        func = @(x,g) dpa.mergeResiduals(Opt,fun,resDpa,x,g);
        x0 = [Path.varAll(:,i);Path.lAll(i)];
        dsBif = mean(diff(Path.sAll(i+(-1:1))));
        OptTrace.dscale0 = max(abs([x0;Opt.g0]),10^-8*ones(numel(x0)+1,1));
        [vari,li,~,~,si] = continuation(func,x0,Opt.g0,Opt.gTarget,dsBif,'Opt',OptTrace);
        PathBifs.varAll = [PathBifs.varAll,NaN(Info.nv+1,1),vari];
        PathBifs.lAll = [PathBifs.lAll,NaN,li];
        PathBifs.sAll = [PathBifs.sAll,NaN,si];
    end
    %% output:
    PathTemp = Path;
    if ~isempty(PathBifs.varAll)
        Path.varAll = [PathTemp.varAll,PathBifs.varAll(1:Info.nv,:)];
        Path.lAll = [PathTemp.lAll,PathBifs.varAll(Info.nv+1,:);
                      Opt.g0*ones(size(PathTemp.lAll)),PathBifs.lAll];
        Path.sAll = [PathTemp.sAll,PathBifs.sAll];
    else
        Path.lAll = [PathTemp.lAll;Opt.g0*ones(size(PathTemp.lAll))];
    end
end