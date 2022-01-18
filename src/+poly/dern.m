%% path continuation - poly.dern
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.10.2020 - Alwin Förster
%
function dp = dern(p,nd)
    try
        no = (length(p)/nd)-1;
        P = reshape(p,[no+1,nd])';
    catch
        error('nd does not match size of p');
    end
    if no>=1
        dP = P(:,1:no).*kron(ones(nd,1),no:-1:1);
    else
        dP = zeros(nd,1);
    end
    dp = reshape(dP',[1,nd*max([1,no])]);
end