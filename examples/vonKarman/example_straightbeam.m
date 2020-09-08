%% Mode selection on Von Karman Straight Beam example
% This example is described in Section 6.1 (Figure 5) of the following 
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
epsilon = 1e-3; % beam thickness-to-length ratio 
midp_height = 0; % height of midpoint relative to the ends (measure of curvature of the beam)
[Geo] = Geometry(epsilon,midp_height);

nElements = 10;
BC = 'B'; % boundary condition: simply supported

%% build finite element model
[Misc,model] = Beam_Model(Geo,nElements,BC); 

n = length(model.freeDOFs); % total number of degrees of freedom
plot_dof = 9; % the degree of freedom at which the response is plotted
%% Mass Matrix
M = model.M(model.freeDOFs,model.freeDOFs);

%% Stiffness Matrix
K = model.K(model.freeDOFs,model.freeDOFs);

%% Damping Matrix
C = model.C(model.freeDOFs,model.freeDOFs);

%% Steady-State tool
SS = SSR(M,C,K,1:n);   % full system

%% Mode selection

% Input parameters for mode selection
param.curvtolerance = 0.05;
param.nmodes = 10;
param.type = 0; % set to 0 if param.nmodes is the max number of modes,
                % set to 1 if param.nmodes is the desired number of modes
param.example = 'simple'; % always use simple if enough memory is available to store X below
param.initialmodes = [1 2 3 4 5]; % optional input modes
X = S11(model); % this conains the quadratic nonlinearities for all dofs (X_i).
                % X_i is symmetric, and can be found from the form M\ddot{q} + C\dot{q} +
                % Kq +\Gamma(q) = F(t) (i denotes the specific row of this equation)

% master mode sets considered
I_1 = 1:10; 
I_2 = modeselect(SS,param,X);
SS1 = SSR(M,C,K,I_1);  % reduced to I_1   
SS2 = SSR(M,C,K,I_2);  % reduced to I_2

%% Incrementing amplitude
T = 2*pi/26;
A = [10:5:35 40:4:172 174 176:1:200 200:1:230]* 10^(-2); % amplitude vector
SS.S = @(x) NonlinearityVK(model,x);
SS.DS = @(x) NonlinearityJacobianVK(model,x); 

tol = 1e-6; % Newton Raphson error tolerance
maxiter = 60;
N2 = nan(length(A),1);
N1 = nan(length(A),1);
NF = nan(length(A),1);
NL = nan(length(A),1);

for j = 1:length(A)
    Fext = @(t,T) A(j)*model.loads(model.freeDOFs) * sin(2*pi*t/T);
    SS.f = Fext;          
    SS.T = T;  
    SS.n_steps = 100;
    SS1 = copyfromfull(SS1,SS);
    SS2 = copyfromfull(SS2,SS);

    if j == 1
        [x2,xd2] = SS1.NewtonRaphson('tol',tol,'maxiter',maxiter); % find out response using Picard iteration
    else
        [x2,xd2] = SS1.NewtonRaphson('init',x2,'tol',tol,'maxiter',maxiter);
    end
    N1(j) = SS1.dt*norm([x2; xd2],'fro');

    if j == 1
        [x3,xd3] = SS2.NewtonRaphson('tol',tol,'maxiter',maxiter); % find out response using Picard iteration
    else
        [x3,xd3] = SS2.NewtonRaphson('init',x3,'tol',tol,'maxiter',maxiter);
    end
    N2(j) = SS2.dt*norm([x3; xd3],'fro');


    if j == 1
        [xf,xdf] = SS.NewtonRaphson('tol',tol,'maxiter',maxiter);
    else
        [xf,xdf] = SS.NewtonRaphson('init',xf,'tol',tol,'maxiter',maxiter);
    end

    NF(j) = SS.dt*norm([xf; xdf],'fro');

    [x_lin, xd_lin] = SS.LinearResponse();
    NL(j) = SS.dt*norm([x_lin; xd_lin],'fro');
    
    disp(vpa(A(j)))
end

figure(3);

%% final amplitude time history
plot(SS.t, xf(plot_dof,:), 'k', 'DisplayName', 'Full','linewidth',1);
axis tight;
legend('show')
ylim([-3.1*10^(-3) 3.1*10^(-3)])
xlabel('time (seconds)');
ylabel('displacement $$w_{4}$$ [m]')
hold on;
plot(SS.t, x_lin(plot_dof,:), 'b', 'DisplayName', 'Linear','linewidth',1);
plot(SS1.t, x2(plot_dof,:),'--g', 'DisplayName', 'Reduced $$I_1$$','linewidth',2);
plot(SS2.t, x3(plot_dof,:),'--r', 'DisplayName', 'Reduced $$I_2$$','linewidth',2);

%% forcing amp - response amp curve
figure
plot(A,NF,'-k','DisplayName', 'Full','linewidth',1); axis tight; hold on;
xlabel('$$F$$ [N]'); ylabel('$$||q||_2$$')
lgd = legend('show');
lgd.Location = 'southeast';
plot(A,NL,'-b', 'DisplayName', 'Linear','linewidth',1);
plot(A,N1,'--g', 'DisplayName', 'Reduced $$I_1$$','linewidth',2);
plot(A,N2,'--r','DisplayName', 'Reduced $$I_2$$','linewidth',2);
ylim([0 0.032])



