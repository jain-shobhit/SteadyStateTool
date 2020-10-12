function Df = SP_nonlinearity_derv_O1(x)
n = length(x)/2;
x = [0; x(n+1:end); 0];
g = 3*(x(2:end) - x(1:end-1)).^2;

D = [-g(2:end-1); 0]; 
E = g(1:end-1)+g(2:end);
F = [0; -g(2:end-1)]; 

Df = spdiags([D E F], -1:1, n,n);
