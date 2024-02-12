clear;
% close all;
clc;
addpath('../../MechConnectX/src/');
addpath('../../HBM/src/');
addpath('../src/');
%% MechConnectX/HBM test case
% MechConnectX
M = 1.234;
K = 170.135;
C = 0.1;
lam1Para = MCXParameter(1);
lam2Para = MCXParameter(1);
lam3Para = MCXParameter(1);
nlE1 = NonlinearElement.Duffing(lam1Para);
nlE2 = NonlinearElement.PST201(lam2Para);
nlE3 = NonlinearElement.PST120(lam3Para);
fNl = NonlinearForce(1);
fNl.addElement(nlE1,1,1);
fNl.addElement(nlE2,1,1);
fNl.addElement(nlE3,1,1);
OmPara = MCXParameter(0);
fEx = HarmonicForce(1,1,[0;0.54673;0],OmPara);
nlS = NonlinearMechanicalSystem(M,K,'damping',C,'nonlinearForce',fNl,...
    'externalForce',fEx);
% HBM
nH = 5;
R = @(Q,Om) resHBM(nlS,nH,Q,Om,OmPara);
%% Continuation
Om0 = 0;
OmE = nlS.eigenAngularFrequencies(end)*2;
Q0 = getQ0(nlS,nH,Om0,OmPara);
ds0 = (OmE-Om0)*10^-3;
dsMax = 1*10^-1;
[Qs,Oms,exitflag,Bifurcation,sAll,Js,breakFunOut,InfoOut] =...
    continuation(R,Q0,Om0,OmE,ds0,...
    'dsMax',dsMax,'jacobianOut','full',...
    'plot','on','checkJacobian','on');
%% post processing
qhs = getAmplitude(nlS,nH,Qs,Oms,OmPara);
sts = hillStability(nlS,nH,Qs,Oms,OmPara,'jacobian',Js,'OmParameter',OmPara);
idxInstab = getInstableSegments(sts);
%% plot
figure(362);
clf;
plot(Oms,qhs,'b-','LineWidth',2);
hold on;
for ii=1:numel(idxInstab(:,1))
    idxii = idxInstab(ii,1):idxInstab(ii,2);
    plot(Oms(idxii),qhs(:,idxii),'r-','LineWidth',2);
end
plot(Oms(Bifurcation.bif(1,:)),qhs(:,Bifurcation.bif(1,:)),'ko','LineWidth',2);
hold off;
try
    xlim([-1,1]*max(abs(Oms(Bifurcation.bif(1,:))-nlS.eigenAngularFrequencies))*1.5+nlS.eigenAngularFrequencies);
catch
    xlim([-1,1]*nlS.eigenAngularFrequencies*0.1+nlS.eigenAngularFrequencies);
end
xticks([0,nlS.eigenAngularFrequencies']);
xticklabels({'0','\omega_0'});
grid on;