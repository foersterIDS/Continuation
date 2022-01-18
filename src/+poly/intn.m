%% path continuation - poly.dern
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.11.2020 - Alwin Förster
%
function ip = intn(p,nd,c)
    try
        no = (length(p)/nd)-1;
        P = reshape(p,[no+1,nd])';
    catch
        error('nd does not match size of p');
    end
    if nargin<=2
        c = zeros(nd,1);
    else
        if numel(c)~=nd
            error('vector of constants must have size [nd x 1]');
        end
        c = c(:);
    end
    if no>=1
        iP = [P.*kron(ones(nd,1),1./((no:-1:0)+1)),c];
    else
        iP = zeros(nd,1);
    end
    ip = reshape(iP',[1,nd*max([1,(no+2)])]);
end