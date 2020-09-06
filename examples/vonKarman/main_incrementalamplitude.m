clear 
clc
close all

%% Von Karman Beam example
epsilon = 1e-3; % h/L 7e-3
midp_height = 0; %0.005 height of midpoint relative to the ends (measure of curvature of the beam)
[Geo] = Geometry(epsilon,midp_height);
nElements = 10;
BC = 'B'; % simply supported
[Misc,model] = Beam_Model(Geo,nElements,BC);

% master mode sets considered
modesel1 = 1:10; 
modesel2 = [1:5 21:22];

n = length(model.freeDOFs);
plot_dof = 9;

%% Mass Matrix
M = model.M(model.freeDOFs,model.freeDOFs);

%% Stiffness Matrix
K = model.K(model.freeDOFs,model.freeDOFs);

%% Damping Matrix
C = model.C(model.freeDOFs,model.freeDOFs);

%% Steady-State tool
SSS = SSR(M,C,K,1:n);
SS1 = SSR(M,C,K,modesel1);   
SS2 = SSR(M,C,K,modesel2);    

T = 2*pi/26;
A = [10:5:35 40:4:172 174 176:1:200 200:1:230]* 10^(-2); % amplitude vector
SSS.S = @(x) NonlinearityVK(model,x);
SSS.DS = @(x) NonlinearityJacobianVK(model,x); 

tol = 1e-6;
maxiter = 60;
N2 = nan(length(A),1);
N1 = nan(length(A),1);
NF = nan(length(A),1);
NL = nan(length(A),1);

for j = 1:length(A)
    Fext = @(t,T) A(j)*model.loads(model.freeDOFs) * sin(2*pi*t/T);
    SSS.f = Fext;          
    SSS.T = T;  
    SSS.n_steps = 100;
    SS1 = copyfromfull(SS1,SSS);
    SS2 = copyfromfull(SS2,SSS);

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
        [xf,xdf] = SSS.NewtonRaphson('tol',tol,'maxiter',maxiter);
    else
        [xf,xdf] = SSS.NewtonRaphson('init',xf,'tol',tol,'maxiter',maxiter);
    end

    NF(j) = SSS.dt*norm([xf; xdf],'fro');

    [x_lin, xd_lin] = SSS.LinearResponse();
    NL(j) = SSS.dt*norm([x_lin; xd_lin],'fro');
    
    disp(vpa(A(j)))
end

figure(2);

%% final amplitude time history
plot(SSS.t, xf(plot_dof,:), 'k', 'DisplayName', 'Full','linewidth',1);
axis tight;
legend('show')
ylim([-3.1*10^(-3) 3.1*10^(-3)])
xlabel('time (seconds)');
ylabel('displacement $$w_{4}$$ [m]')
hold on;
plot(SSS.t, x_lin(plot_dof,:), 'b', 'DisplayName', 'Linear','linewidth',1);
plot(SS1.t, x2(plot_dof,:),'--g', 'DisplayName', 'Reduced $$I_1$$','linewidth',2);
plot(SS2.t, x3(plot_dof,:),'--r', 'DisplayName', 'Reduced $$I_2$$','linewidth',2);

%% forcing amp - response amp curve
% plot(A,NF,'-k','DisplayName', 'Full','linewidth',1); axis tight; hold on;
% xlabel('$$F$$ [N]'); ylabel('$$||q||_2$$')
% lgd = legend('show');
% lgd.Location = 'southeast';
% plot(A,NL,'-b', 'DisplayName', 'Linear','linewidth',1);
% plot(A,N1,'--g', 'DisplayName', 'Reduced $$I_1$$','linewidth',2);
% plot(A,N2,'--r','DisplayName', 'Reduced $$I_2$$','linewidth',2);
% ylim([0 0.032])



