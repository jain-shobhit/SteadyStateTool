%% Mode Selection on Von Karman curved Beam Example
% This example is described in Section 6.2 (Figure 8) of the following 
% article. 
% G. Buza, S. Jain, G. Haller, Using Spectral Submanifolds for Optimal 
% Mode Selection in Model Reduction, (2020) Preprint available on arXiv.org

% The finite element model is taken from the following article. 
% Jain, S., Tiso, P., & Haller, G. (2018). Exact nonlinear model reduction 
% for a von Kármán beam: slow-fast decomposition and spectral submanifolds. 
% Journal of Sound and Vibration, 423, 195–211. 
% https://doi.org/10.1016/J.JSV.2018.01.049

clear 
clc
close all

%% parameters
epsilon = 7e-3; % beam thickness-to-length ratio 
midp_height = 0.005; % height of midpoint relative to the ends (measure of curvature of the beam)
[Geo] = Geometry(epsilon,midp_height);

nElements = 10;
BC = 'B'; % boundary condition: simply supported

%% build finite element model
[Misc,model] = Beam_Model(Geo,nElements,BC); 

n = length(model.freeDOFs); % total number of degrees of freedom

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

%% Steady-State tool
SSfull = SSR(M,C,K,1:n);        % full system

%% Mode selection

% Input parameters for mode selection
param.curvtolerance = 0.05;
param.nmodes = 10;
param.type = 1; % set to 0 if param.nmodes is the max number of modes,
% set to 1 if param.nmodes is the desired number of modes
param.example = 'simple'; % always use simple if enough memory is available to store X below
X = S11(model); % this conains the quadratic nonlinearities for all dofs (X_i).
% X_i is symmetric, and can be found from the form M\ddot{q} + C\dot{q} +
% Kq +\Gamma(q) = F(t) (i denotes the specific row of this equation)

% master mode sets considered
modesel1 = 1:10; 
modesel2 = modeselect(SSfull,param,X);
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



