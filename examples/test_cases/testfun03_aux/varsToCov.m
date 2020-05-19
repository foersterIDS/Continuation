function [Kxx,mux] = varsToCov(nlsx,vars)
    mf = nlsx.isMeanFree;
    nx = nlsx.nx;
    Cov_info = nlsx.fnl.Cov_info;
    nK = length(Cov_info);
    Kxx = sparse(Cov_info-floor(Cov_info/nx-10^-8)*nx,ceil(Cov_info/nx),vars(1:nK),nx,nx);
    Kxx = Kxx'+Kxx-diag(diag(Kxx));
    if mf
        mux = zeros(nx,1);
    else
        Mux_info = nlsx.fnl.Mux_info;
        nM = length(Mux_info);
        mux = sparse(Mux_info,1,vars(nK+(1:nM)),nx,1);
    end
end