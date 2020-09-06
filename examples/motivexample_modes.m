clear 
clc
close all

%% Mass-spring example
M = eye(3);
D1 = 0.01;
D2 = 0.2;
D3 = 0.08; 
ome1 = 2;
ome2 = 3;
ome3 = 5;
K = diag([ome1^2 ome2^2 ome3^2]);
C = 2*diag([ome1*D1 ome2*D2 ome3*D3]);

% Nonlinearity
S = @(x)SP_nonlinearity3(x,ome1,ome2,ome3);
DS = @(x)SP_nonlinearity_derv3(x,ome1,ome2,ome3);

% Loading
A = 0.02;                       % loading amplitude       
f0 = [1; 0; 0];ones(3,1);                 % loading shape
f = @(t,T)A*f0*cos(2*pi*t/T);   % loading function

%% SSR Package
SS = SSR(M,C,K,1);

% SSM computation
R = cell(1,2);
R{1} = [ome2^2/2];
R{2} = [ome3^2/2];
W = SS.SSM2(R);

% Display mode recommendation
disp(norm(W{1}))
disp(norm(W{2}))

function f = SP_nonlinearity3(x,o1,o2,o3)
f = [o1^2/2*(3*x(1)^2+x(2)^2+x(3)^2)+o2^2*x(1)*x(2)+o3^2*x(1)*x(3)+(o1^2+o2^2+o3^2)/2*x(1)*(x(1)^2+x(2)^2+x(3)^2);...
    o2^2/2*(3*x(2)^2+x(1)^2+x(3)^2)+o1^2*x(1)*x(2)+o3^2*x(2)*x(3)+(o1^2+o2^2+o3^2)/2*x(2)*(x(1)^2+x(2)^2+x(3)^2);...
    o3^2/2*(3*x(3)^2+x(2)^2+x(1)^2)+o2^2*x(3)*x(2)+o1^2*x(1)*x(3)+(o1^2+o2^2+o3^2)/2*x(3)*(x(1)^2+x(2)^2+x(3)^2) ];
end

function Df = SP_nonlinearity_derv3(x,o1,o2,o3)
Df = zeros(3);
Df(1,1) = 3*o1^2*x(1)+o2^2*x(2)+o3^2*x(3) + (o1^2+o2^2+o3^2)/2*(3*x(1)^2 + x(2)^2+x(3)^2);
Df(2,2) = 3*o2^2*x(2)+o1^2*x(1)+o3^2*x(3) + (o1^2+o2^2+o3^2)/2*(3*x(2)^2 + x(1)^2+x(3)^2);
Df(3,3) = (3*o3^2*x(3)+o2^2*x(2)+o1^2*x(1) + (o1^2+o2^2+o3^2)/2*(3*x(3)^2 + x(2)^2+x(1)^2));
Df(1,2) = o1^2*x(2)+o2^2*x(1)+(o1^2+o2^2+o3^2)*x(1)*x(2);
Df(1,3) = o1^2*x(3)+o3^2*x(1)+(o1^2+o2^2+o3^2)*x(1)*x(3);
Df(2,3) = o2^2*x(3)+o3^2*x(2)+(o1^2+o2^2+o3^2)*x(2)*x(3);
Df(3,2) = Df(2,3);
Df(2,1) = Df(1,2);
Df(3,1) = Df(1,3);
end