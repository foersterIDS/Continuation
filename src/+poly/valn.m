%% path continuation - poly.valn
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   12.10.2020 - Alwin Förster
%
function yp = valn(p,x,nd)
    x = x(:).';
    nx = length(x);
    try
        no = (length(p)/nd)-1;
        P = reshape(p,[no+1,nd])';
    catch
        error('nd does not match size of p');
    end
    yp = zeros(nd,nx);
    for i=0:no
        yp = yp+P(:,i+1)*x.^(no-i);
    end
end