function f = SP_nonlinearity_O1(x)
n = numel(x)/2;
x = [0; x(n+1:end); 0];
g = (x(2:end) - x(1:end-1)).^3;

f = g(1:end-1)-g(2:end);
