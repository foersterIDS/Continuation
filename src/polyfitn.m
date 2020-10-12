%% path continuation - polyfitn
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   12.10.2020 - Alwin Förster
%
function p = polyfitn(x,y,no)
    [nd,nx] = size(y);
    if nx<(no+1)
        error('polynomial order to high');
    end
    XX = NaN(nx,no+1);
    for i=0:no
        XX(:,i+1) = x(:).^i;
    end
    X = kron(eye(nd),XX);
    yy = NaN(nx*nd,1);
    for i=1:nd
        yy(nx*(i-1)+(1:nx)) = y(i,:).';
    end
    p = ((X'*X)\X'*yy).';
end