function [K] = kette(k,Ndof,fixed)
    if nargin<3
        fixed = 1; % DOf mit Umgebungkopplung
    end
    K = zeros(Ndof);
    if Ndof==1
        K = k;
    else
        for i=1:Ndof
            if i==1
                K(1,1) = k;
                K(1,2) = -k;
            elseif i==Ndof
                K(Ndof,Ndof-1) = -k;
                K(Ndof,Ndof) = k;
            else
                K(i,i-1) = -k;
                K(i,i) = 2*k;
                K(i,i+1) = -k;
            end
            if sum(i==fixed)
                K(i,i) = K(i,i)+k;
            end
        end
    end
end