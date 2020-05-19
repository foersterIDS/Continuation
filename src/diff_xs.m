function [dxs] = diff_xs(xs)
    if length(xs(1,:))==1
        n = length(xs(:,1));
        dxs = [zeros(n-1,1);1];
    elseif length(xs(1,:))>1
        dxs = xs(:,end)-xs(:,end-1);
    else
        error('xs is empty');
    end
end