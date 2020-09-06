clear 
clc
close all

%% Von Karman Beam example
epsilon = 7e-3; % h/L 7e-3 1e-3
midp_height = 0.005; %0.005 height of midpoint relative to the ends (measure of curvature of the beam)
[Geo] = Geometry(epsilon,midp_height);
nElements = 10;
BC = 'B'; % simply supported
[Misc,model] = Beam_Model(Geo,nElements,BC);

modesel1 = [1:10];
modesel2 = [1:8 12 17];
n = length(model.freeDOFs);

%% Mass Matrix
M = model.M(model.freeDOFs,model.freeDOFs);

%% Stiffness Matrix
K = model.K(model.freeDOFs,model.freeDOFs);

%% Damping Matrix
C = model.C(model.freeDOFs,model.freeDOFs);

%% Nonlinearity
S =@(x) NonlinearityVK(model,x);
DS =@(x) NonlinearityJacobianVK(model,x);

A = 80; % amplitude for which the FRC is computed
f = @(t,T) A*model.loads(model.freeDOFs) * sin(2*pi*t/T); % loading function

%% SSR package
SSfull = SSR(M,C,K,1:n);        % full system
SSred1 = SSR(M,C,K,modesel1);   % reduced to I_1
SSred2 = SSR(M,C,K,modesel2);   % reduced to I_2

SSfull.S = S;
SSfull.DS = DS;
SSfull.f = f;
SSred1 = copyfromfull(SSred1,SSfull);
SSred2 = copyfromfull(SSred2,SSfull);

%% Sequential continuation
omega1 = SSfull.omega(1);
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



