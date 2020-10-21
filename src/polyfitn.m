%% path continuation - polyfitn
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   12.10.2020 - Alwin Förster
%
function p = polyfitn(x,y,no)
    warning('');
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
    p = ((X'*X)\X'*yy).';
    [warnMsg, warnId] = lastwarn;
    if ~isempty(warnMsg) && no>1
        p = polyfitn(x((end-no+1):end),y(:,(end-no+1):end),no-1);
    end
end