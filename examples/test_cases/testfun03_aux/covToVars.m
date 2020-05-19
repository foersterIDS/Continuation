function [vars] = covToVars(nlsx,Kxx,mux)
    mf = nlsx.isMeanFree;
    if mf
        Cov_info = nlsx.fnl.Cov_info;
        vars = Kxx(Cov_info);
    else
        Cov_info = nlsx.fnl.Cov_info;
        Mux_info = nlsx.fnl.Mux_info;
        vars = [Kxx(Cov_info);mux(Mux_info)];
    end
end