%% path continuation - dpa.trace
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   04.04.2022 - Alwin Förster
%
function trace(fun,dpaPoints,Info,Opt,Path)
    %% settings for dpa:
    dpaTrace = dpaPoints;
    npoints = numel(dpaTrace);
    OptTrace = Opt;
    OptTrace.direction = [zeros(size(OptTrace.direction));sign(Opt.gTarget-Opt.g0)];
    OptTrace.l0 = Opt.g0;
    OptTrace.lTarget = Opt.gTarget;
    OptTrace.dpaGammaVar = true;
    PathDpa = struct('varAll',[],'lAll',[],'sAll',[]);
    %% start dpa:
    for ii=1:npoints
        i = dpaTrace(ii);
        func = @(x,g) dpa.mergeResiduals(Opt,fun,Opt.dpaResidual,x,g);
        x0 = [Path.varAll(:,i);Path.lAll(i)];
        dsDpa = max([mean(diff(Path.sAll(i+(-1:1))))/100,Opt.dsMin]);
        OptTrace.dscale0 = max(abs([x0;Opt.g0]),10^-8*ones(numel(x0)+1,1));
        [varI,lI,~,~,sI] = continuation(func,x0,Opt.g0,Opt.gTarget,dsDpa,'Opt',OptTrace);
        PathDpa.varAll = [PathDpa.varAll,NaN(Info.nv+1,1),varI];
        PathDpa.lAll = [PathDpa.lAll,NaN,lI];
        PathDpa.sAll = [PathDpa.sAll,NaN,sI];
    end
    %% output:
    PathTemp = Path;
    if ~isempty(PathDpa.varAll)
        Path.varAll = [PathTemp.varAll,PathDpa.varAll(1:Info.nv,:)];
        Path.lAll = [PathTemp.lAll,PathDpa.varAll(Info.nv+1,:);
                      Opt.g0*ones(size(PathTemp.lAll)),PathDpa.lAll];
        Path.sAll = [PathTemp.sAll,PathDpa.sAll];
    else
        Path.lAll = [PathTemp.lAll;Opt.g0*ones(size(PathTemp.lAll))];
    end
end