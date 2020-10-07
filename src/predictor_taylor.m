%% path continuation - predictor_taylor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   05.10.2020 - Alwin Förster
%
function [xp,ch] = predictor_taylor(vars,ls,nd,ds)
    nd = min([length(ls)-1,nd]);
    xs = [vars;ls];
    nx = length(xs(:,end));
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
    
    if length(ls) > 2
        vh = height(vars);
        theta = 0;
        
        % calculate greatest angle
        for k = 1:vh
            vr = vars(k,end-2:end);
            lr = ls(end-2:end);
            vec = [lr;vr];
            v1 = vec(:,end) - vec(:,end-1);
            v2 = vec(:,end-1) - vec(:,end-2);
            
            thetn = real(acos(max(min(dot(v1,v2)/(norm(v1)*norm(v2)),1),-1)));

            if thetn > theta
                theta = thetn;
            end
        end
        
        % define highest value for angle
        theta_set = 5 * pi/180;
        
        % calculate ratio
        ch = theta / theta_set;
        
        if ch < 1 % if ratio < 1, set ch to 1
            ch = 1;
        elseif ch > 1.5 % high ratios are punished stronger
            ch = 5;
        end
    else
        ch = 1;
    end

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