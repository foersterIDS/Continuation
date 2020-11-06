clear all;
clc;

x0 = [500;10^-2];
sc1 = [1;1]; % keine Skalierung
sc2 = [500;10^-2]; % Skalierung in der Groessenordnung
fun1 = @(x) res_fun(x,sc1);
fun2 = @(x) res_fun(x,sc2);
fun3 = @(x) res_fun(x,ones(2,1));
x0sc1 = x0./sc1;
x0sc2 = x0./sc2;
x03 = x0;

opt12 = optimoptions('fsolve','SpecifyObjectiveGradient',false); % test analytic/numeric jacobian
[xsc1,~,exitflag1,output1,Jsc1] = fsolve(fun1,x0sc1,opt12);
x1 = xsc1.*sc1;
[xsc2,~,exitflag2,output2,Jsc2] = fsolve(fun2,x0sc2,opt12);
x2 = xsc2.*sc2;

opt3 = optimoptions('fsolve','SpecifyObjectiveGradient',true,'TypicalX',sc2); % test analytic/numeric jacobian
[x3,~,exitflag3,output3,J3] = fsolve(fun3,x03,opt3);

[~,Janasc1] = fun1(xsc1);
[~,Janasc2] = fun2(xsc2);

clc;
fprintf('-- Skalierung #1 --\n');
fprintf('exitflag: %d\n',exitflag1);
fprintf('x1 = %.2e\t|\tx2 = %.2e\n',x1(1),x1(2));
fprintf('iterations: %d\n',output1.iterations);
fprintf('jac-error: %.2e\n',norm(Janasc1(:)-Jsc1(:))/norm(Janasc1(:)));
fprintf('\n');
fprintf('-- Skalierung #2 --\n');
fprintf('exitflag: %d\n',exitflag2);
fprintf('x1 = %.2e\t|\tx2 = %.2e\n',x2(1),x2(2));
fprintf('iterations: %d\n',output2.iterations);
fprintf('jac-error: %.2e\n',norm(Janasc2(:)-Jsc2(:))/norm(Janasc2(:)));
fprintf('\n');
fprintf('-- Skalierung #3 --\n');
fprintf('exitflag: %d\n',exitflag3);
fprintf('x1 = %.2e\t|\tx2 = %.2e\n',x3(1),x3(2));
fprintf('iterations: %d\n',output3.iterations);
fprintf('jac-error: %.2e\n',norm(Janasc1(:)-J3(:))/norm(Janasc1(:)));