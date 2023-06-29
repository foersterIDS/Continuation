%% path continuation - poly.shiftn
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.06.2023 - Alwin Förster
%
function ps = shiftn(p,xs,nd)
    %% arguments
    arguments
        p (1,:) double
        xs (1,1) double
        nd (1,1) double
    end
    %% reshape (p --> P)
    try
        no = (length(p)/nd)-1;
        P = reshape(p,[no+1,nd]).';
    catch
        error('nd does not match size of p');
    end
    %% shift
    Ps = zeros(nd,no+1);
    for ii=0:no
        for jj=0:ii
            kk = (no+1)-(ii-jj);
            bc = poly.binomialCoefficient(ii,jj);
            Ps(:,kk) = Ps(:,kk)+bc*P(:,no+1-ii)*xs^jj;
        end
    end
    %% reshape (Ps --> ps)
    ps = reshape(Ps.',[1,nd*(no+1)]);
end