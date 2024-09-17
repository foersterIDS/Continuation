%% path continuation - bifurcation.parameterTrace
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   31.03.2022 - Alwin Förster
%
function parameterTrace(oih,fun)
    %% settings for dpa:
    bifTrace = oih.bifurcation.bif(1,:);
    nbif = numel(bifTrace);
    optTrace = oih.opt;
    optTrace = aux.setoff(optTrace,'bifurcation');
    fnOiS = fieldnames(oih.optIsSet);
    optTrace = rmfield(optTrace,fnOiS(~cell2mat(struct2cell(oih.optIsSet)))); 
    optTrace.direction = [zeros(size(oih.opt.direction));sign(oih.opt.gTarget-oih.opt.g0)];
    optTrace.l0 = oih.opt.g0;
    optTrace.lTarget = oih.opt.gTarget;
    optTrace.dpaGammaVar = true;
    varAllBifs = [];
    lAllBifs = [];
    %% start dpa:
    for ii=1:nbif
        idx = bifTrace(ii);
        sc = oih.bifurcation.scaling(ii);
        resDpa = @(v,l,g) dpa.resBif(fun,[v;l],g,oih,sc);
        func = @(x,g) dpa.mergeResiduals(oih,fun,resDpa,x,g);
        x0 = [oih.path.varAll(:,idx);oih.path.lAll(idx)];
        dsBif = mean(diff(oih.path.sAll(idx+(-1:1))));
        optTrace.dscale0 = max(abs([x0;oih.opt.g0]),10^-8*ones(numel(x0)+1,1));
        [vari,li,~,~,si] = continuation(func,x0,oih.opt.g0,oih.opt.gTarget,dsBif,'Opt',optTrace);
        varAllBifs = [varAllBifs,NaN(oih.info.nv+1,1),vari];
        lAllBifs = [lAllBifs,NaN,li];
    end
    %% output:
    if ~isempty(varAllBifs)
        varAllTemp = [oih.path.varAll,varAllBifs(1:oih.info.nv,:)];
        lAllTemp = [oih.path.lAll,varAllBifs(oih.info.nv+1,:);
                    oih.opt.g0*ones(size(oih.path.lAll)),lAllBifs];
        oih.path.overwrite(varAllTemp,lAllTemp,[],true);
    else
        lAllTemp = [oih.path.lAll;oih.opt.g0*ones(size(oih.path.lAll))];
        oih.path.overwrite([],lAllTemp,[],false);
    end
end