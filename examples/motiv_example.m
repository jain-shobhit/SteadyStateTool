%% Mass-spring motivational example for mode selection

% This example is described in Section 4.3 (Figure 2) of the following 
% article 
% G. Buza, S. Jain, G. Haller, Using Spectral Submanifolds for Optimal 
% Mode Selection in Model Reduction, (2020) Preprint available on arXiv.org

clear 
clc
close all

%% properties
M = eye(3);
D1 = 0.01;
D2 = 0.02;
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
A = 0.02;                     % loading amplitude       
f0 = [1;0;0];                 % loading shape
f = @(t,T)A*f0*sin(2*pi*t/T); % loading function

%% Mode selection

SS = SSR(M,C,K,1);

% SSM computation
R = cell(1,2);
R{1} = ome2^2/2;
R{2} = ome3^2/2;
W = SS.SSM2(R);
disp(['norm(W_2)=' num2str(norm(W{1})) ', norm(W_3)=' num2str(norm(W{2}))])
[~, ind] = max([norm(W{1}) norm(W{2})]);

% master mode sets considered
I_1 = [1 2];        % mode selection 1
I_2 = [1 ind+1];    % mode selection 2

%% SSR Package
SSfull = SSR(M,C,K,[1 2 3]);   % full system
SSred1 = SSR(M,C,K,I_1);     % reduced to I_1
SSred2 = SSR(M,C,K,I_2);     % reduced to I_2
SSfull.S = S;
SSfull.DS = DS; 
SSfull.f = f;
SSred1 = copyfromfull(SSred1,SSfull);
SSred2 = copyfromfull(SSred2,SSfull);

%% Sequential continuation
omega1 = min(SSfull.omega);
Omega_range = [0.7 1.3]*omega1;
T_range = 2*pi./Omega_range;
ncontsteps = 100;

[Omega_arrayf,Nf] = SSfull.FRC(T_range,ncontsteps);
[Omega_array1,N1] = SSred1.FRC(T_range,ncontsteps);
[Omega_array2,N2] = SSred2.FRC(T_range,ncontsteps);

Nlin = nan(length(Omega_arrayf),1);
for j = 1:length(Omega_arrayf)
    SSfull.T = 2*pi/Omega_arrayf(j);
    [x_lin, xd_lin] = SSfull.LinearResponse();
    LinSol = [x_lin; xd_lin];
    Nlin(j) = SSfull.dt*norm(LinSol,'fro');
end

figure(2); semilogy(Omega_arrayf,Nf,'-k','DisplayName', 'Full','linewidth',1); axis tight; grid on; hold on;
xlabel('$$\Omega$$ [rad/s]'); ylabel('$$||q||_2$$')
legend('show')
semilogy(Omega_arrayf,Nlin,'-b', 'DisplayName', 'Linear','linewidth',1);
semilogy(Omega_array1,N1,'--g', 'DisplayName', 'Reduced $$I_1$$','linewidth',2);
semilogy(Omega_array2,N2,'--r','DisplayName', 'Reduced $$I_2$$','linewidth',2);


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