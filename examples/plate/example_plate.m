%% Steady state response of the curved plate example
% This example is described in Section 5.3 (Figure 5) of the following 
% article. 
% G. Buza, S. Jain, G. Haller, Integral Equations & Model Reduction For 
% Fast Computation of Nonlinear Periodic Response, (2020) Preprint available on arXiv.org

% The finite element model is taken from the following article.
% Tiso, P. Finite element based reduction methods for static and dynamic 
% analysis of thin-walled structures
% PhD thesis (Delft University of Technology, 2006).

close all
clear
clc

%% parameters
pars = parameters();
[model,Misc,pars] = Initialize(pars);
nDOFs = Misc.ndof;

%% FE matrices
activeDOFs = Misc.freedofs;
M = Misc.Mm(activeDOFs,activeDOFs);
K = Misc.Km(activeDOFs,activeDOFs);
C = Misc.Cm(activeDOFs,activeDOFs);

%% plot undeformed mesh
figure(1)
def_plot(1e-16+zeros(nDOFs,1),0.1,Misc.nodes,Misc.elements)

%% Nonlinearity
S = @(x) Nonlinearity(model,Misc,x);
DS = @(x) NonlinearityJacobian(model,Misc,x);

amplitudes = [1 2 4]*10^(-2); % amplitudes for which FRCs are computed
ncontsteps = 40;              % number of continuation steps

%% Mode selection

% Input parameters for mode selection
modeselpars.curvtolerance = 0.1;
modeselpars.nmodes = 10;
modeselpars.type = 1; % set to 0 if param.nmodes is the max number of modes,
% set to 1 if param.nmodes is the desired number of modes

% use this approach if the coefficients cannot be stored for the full
% system (compare with modeselect.m)
modeselpars.example = 'plate'; 
modeselpars.fulldof = nDOFs;
modeselpars.freedofs = Misc.freedofs;
modeselpars.model = model;

System.M=M;System.C=C;System.K=K;System.sys_order='second';
full = SSR(System);
mode_choice = modeselect(full,modeselpars);

%% Steady-State tool
SSfull = full;
SSred = SSR(System,mode_choice); 

SSred.S = S;
SSred.DS = DS;
SSfull.S = S;
SSfull.DS = DS;

%% Sequential continuation
omega1 = SSfull.omega(2);
Omega_range = [0.8 1.07]*omega1;
T_range = 2*pi./Omega_range;

figure(2)
for i = 1:length(amplitudes)
    f = @(t,T) amplitudes(i)*Misc.loads(activeDOFs) * sin(2*pi*t/T);
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