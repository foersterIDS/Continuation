function [R] = residual_filter2(nlsx,fex,vars,om)
    fex.fi.updateOm(om);
    nlsz = fex.applyExcitation(nlsx);
    nlsz.setMeanFree(nlsx.isMeanFree);
    [Kzzi,muzi] = varsToCov(nlsz,vars);
    lsz = nlsz.getLinearizedSystem(muzi,Kzzi);
    muzip1 = nlsz.getMeanX(Kzzi);
    Kzzip1 = lsz.getKXX();
    varsip1 = covToVars(nlsz,Kzzip1,muzip1);
    R = vars-varsip1;
end