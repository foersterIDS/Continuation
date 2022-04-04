%% path continuation - poly.fitn
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   12.10.2020 - Alwin F�rster
%
function p = fitn(x,y,no,enforce_no)
    [nd,nx] = size(y);
    if nx<(no+1)
        error('polynomial order to high');
    end
    if nargin<=3
        enforce_no = false;
    end
    X = kron(eye(nd),x(:).^(0:no));
    yT = y.';
    yy = yT(:);
    if rcond(X'*X)<10^-15 && no>1 && ~enforce_no
        pr = poly.fitn(x((end-no+1):end),y(:,(end-no+1):end),no-1);
        p = zeros(1,nd*(no+1));
        ind = 1:(nd*(no+1));
        ind((no+1):(no+1):end) = [];
        p(ind) = pr;
    else
        p = ((X'*X)\X'*yy).';
    end
    i = 1:((no+1)*nd);
    j = ceil((1:(nd*(no+1)))/(no+1));
    ind = (j-1)*(no+1)-(no+1)*floor(no*i/(no+1))+no*i+1;
    p = p(ind);
end