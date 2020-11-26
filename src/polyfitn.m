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
    yT = y.';
    yy = yT(:);
    if rcond(X'*X)<10^-15 && no>1
        pr = polyfitn(x((end-no+1):end),y(:,(end-no+1):end),no-1);
        p = zeros(1,nd*(no+1));
        ind = 1:(nd*(no+1));
        ind((no+1):(no+1):end) = [];
        p(ind) = pr;
    else
        p = ((X'*X)\X'*yy).';
    end
    i = 1:((no+1)*nd);
    j = kron(1:nd,ones(1,no+1));
    ind = (j-1)*(no+1)-(no+1)*floor(no*i/(no+1))+no*i+1;
    p = p(ind);
end