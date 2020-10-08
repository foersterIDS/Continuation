%% path continuation - predictor_taylor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   05.10.2020 - Alwin Förster
%
function [xp] = predictor_taylor(vars,ls,nd,ds)
    nd = min([length(ls)-1,nd]);
    xs = [vars;ls];
    nx = length(xs(:,end));
    %% calc arc-length approximation:
    s = zeros(nd,1);
    for i=1:nd
        if i==1
%             s(i) = sqrt(sum((xs(:,end)-xs(:,end-1)).^2));
            s(i) = norm(xs(:,end)-xs(:,end-1));
        else
%             s(i) = s(i-1)+sqrt(sum((xs(:,end-i+1)-xs(:,end-i)).^2));
            s(i) = s(i-1)+norm(xs(:,end-i+1)-xs(:,end-i));
        end
    end
    %% calc taylor-predictor:
    f = kron(ones(nd,1),xs(:,end));
    fdi = zeros(nd*nx,1);
    Mdi = zeros(nd*nx);
    Mpi = zeros(nx,nd*nx);
    for i=1:nd
        fdi(nx*(i-1)+(1:nx),1) = xs(:,end-i);
        Mdi(:,nx*(i-1)+(1:nx)) = kron((-s).^i./factorial(i),eye(nx));
        Mpi(:,nx*(i-1)+(1:nx)) = ds.^i./factorial(i)*eye(nx);
    end
    dfdi = Mdi\(fdi-f);
    xp = xs(:,end)+Mpi*dfdi;
end