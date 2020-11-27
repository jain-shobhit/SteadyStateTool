%% Steady state response of the Von Karman curved Beam Example
% This example is described in Section 5.2 (Figure 3) of the following 
% article. 
% G. Buza, S. Jain, G. Haller, Integral Equations & Model Reduction For 
% Fast Computation of Nonlinear Periodic Response, (2020) Preprint available on arXiv.org

% The finite element model is taken from the following article. 
% Jain, S., Tiso, P., & Haller, G. (2018). Exact nonlinear model reduction 
% for a von K�rm�n beam: slow-fast decomposition and spectral submanifolds. 
% Journal of Sound and Vibration, 423, 195�211. 
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

%% modes selected for a projection based ROM, see example_curvedbeam for details
mode_choice = [1:7 12 17]; 

%% Mass Matrix
M = model.M(model.freeDOFs,model.freeDOFs);

%% Stiffness Matrix
K = model.K(model.freeDOFs,model.freeDOFs);

%% Damping Matrix
C = model.C(model.freeDOFs,model.freeDOFs);

%% Nonlinearity
S = @(x) NonlinearityVK(model,x);
DS = @(x) NonlinearityJacobianVK(model,x);

amplitudes = [1 2 4 8 ]*10^(1); % amplitudes for which FRCs are computed
ncontsteps = 100;               % number of continuation steps

%% Steady-State tool
System.M=M;System.C=C;System.K=K;System.sys_order='second';
SSfull = SSR(System);               % full system
SSred = SSR(System,mode_choice);    % reduced system

SSred.S = S;
SSred.DS = DS;
SSfull.S = S;
SSfull.DS = DS;

%% Sequential continuation
omega1 = SSfull.omega(1);
Omega_range = [0.7 1.3]*omega1;             
T_range = 2*pi./Omega_range;

for i = 1:length(amplitudes)
    f = @(t,T) amplitudes(i)*model.loads(model.freeDOFs) * sin(2*pi*t/T);
    SSred.f = f;
    SSfull.f = f;
    
    [Omega_arrayori,Nori] = SSfull.FRC(T_range,ncontsteps,'original');  
    [Omega_array_r,Nr] = SSfull.FRC(T_range,ncontsteps,'reformulated');
    [Omega_array_rr,Nrr] = SSred.FRC(T_range,ncontsteps,'reformulated');  
    
    semilogy(Omega_arrayori,Nori,'-k','linewidth',2); hold on
    semilogy(Omega_array_r,Nr,'--r','linewidth',2)
    semilogy(Omega_arrayori(1:10:end),Nrr(1:10:end),'sb','MarkerFaceColor','b');    
end

axis tight; grid on;
legend('original full','reformulated full','reformulated reduced','location','SW')
xlabel('excitation frequency $$\Omega$$');
ylabel('response $$||q||_2$$')
yLimits = get(gca,'YLim');
ylim([yLimits(1) 0.03])
xlim(Omega_range)
annotation('textbox', [0.77 0.285 .07 .05],'String',' $$F=10$$ ','LineStyle','none','FontSize',14,'interpreter','latex');
annotation('textbox', [0.77 0.44 .07 .05],'String',' $$F=20$$ ','LineStyle','none','FontSize',14,'interpreter','latex');
annotation('textbox', [0.77 0.6 .07 .05],'String',' $$F=40$$ ','LineStyle','none','FontSize',14,'interpreter','latex');
annotation('textbox', [0.77 0.74 .07 .05],'String',' $$F=80$$ ','LineStyle','none','FontSize',14,'interpreter','latex');
