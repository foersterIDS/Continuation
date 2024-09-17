%% path continuation - dpa.trace
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   04.04.2022 - Alwin Förster
%
function trace(fun,dpaPoints,oih)
    %% settings for dpa:
    dpaTrace = dpaPoints;
    npoints = numel(dpaTrace);
    optTrace = oih.opt;
    fnOiS = fieldnames(oih.optIsSet);
    optTrace = rmfield(optTrace,fnOiS(~cell2mat(struct2cell(oih.optIsSet)))); 
    optTrace.direction = [zeros(size(oih.opt.direction));sign(oih.opt.gTarget-oih.opt.g0)];
    optTrace.l0 = oih.opt.g0;
    optTrace.lTarget = oih.opt.gTarget;
    optTrace.dpaGammaVar = true;
    varAllDpa = [];
    lAllDpa = [];
    sAllDpa = [];
    %% start dpa:
    for ii=1:npoints
        idx = dpaTrace(ii);
        func = @(x,g) dpa.mergeResiduals(oih,fun,oih.opt.dpaResidual,x,g);
        x0 = [oih.path.varAll(:,idx);oih.path.lAll(idx)];
        dsDpa = max([mean(diff(oih.path.sAll(idx+(-1:1))))/100,oih.opt.dsMin]);
        optTrace.dscale0 = max(abs([x0;oih.opt.g0]),10^-8*ones(numel(x0)+1,1));
        [varI,lI,~,~,sI] = continuation(func,x0,oih.opt.g0,oih.opt.gTarget,dsDpa,'Opt',optTrace);
        varAllDpa = [varAllDpa,NaN(oih.info.nv+1,1),varI];
        lAllDpa = [lAllDpa,NaN,lI];
        sAllDpa = [sAllDpa,NaN,sI];
    end
    %% output:
    if ~isempty(varAllDpa)
        varAllTemp = [oih.path.varAll,varAllDpa(1:oih.info.nv,:)];
        lAllTemp = [oih.path.lAll,varAllDpa(oih.info.nv+1,:);
                    oih.opt.g0*ones(size(oih.path.lAll)),lAllDpa];
        oih.path.overwrite(varAllTemp,lAllTemp,[],true);
    elseif ~oih.opt.bifurcation.parameterTrace
        lAllTemp = [oih.path.lAll;oih.opt.g0*ones(1,oih.path.nAll)];
        oih.path.overwrite([],lAllTemp,[],false);
    end
end