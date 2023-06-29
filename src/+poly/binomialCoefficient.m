%% path continuation - poly.binomialCoefficient
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.06.2023 - Alwin Förster
%
function b = binomialCoefficient(n,k)
    if n<0
        n = -n;
        warning('n is negative!');
    end
    if k<0
        k = -k;
        warning('k is negative!');
    end
    fn = factorial(n);
    fk = factorial(k);
    if n>=k
        fnk = factorial(n-k);
        b = fn/(fk*fnk);
    else
        fnk = factorial(k-n);
        b = fk/(fn*fnk);
        warning('n must be larger than k or equal! --> calculated result: k over n');
    end
end