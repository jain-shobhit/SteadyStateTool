function [scal,scaldir] = scalg(A,W)
% This function computes directional scalar curvatures.

m = length(W{1});

SUM = 0;
scaldir = zeros(1,length(W));
for l = 1:length(W)
    Wl = W{l};
    sumk = 0;
    for i = 1:m
        for j = 1:m
            a = 4 * ( Wl(i,i) * Wl(j,j) - (Wl(i,j))^2 +...
                4* (Wl(i,:)*A(:,i))*(Wl(j,:)*A(:,j))-...
                (Wl(i,:)*A(:,j) + Wl(j,:)*A(:,i))^2);
            SUM = SUM + a;
            sumk = sumk+a;
        end
    end
    scaldir(l) = sumk;
end

scal = SUM/(m^2);
end