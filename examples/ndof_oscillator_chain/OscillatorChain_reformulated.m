%% Steady state response of the Von Karman curved Beam Example
% This example is described in Section 5.1 (Figure 2) of the following 
% article. 
% G. Buza, S. Jain, G. Haller, Integral Equations & Model Reduction For 
% Fast Computation of Nonlinear Periodic Response, (2020) Preprint available on arXiv.org

clear
clc
close all

%% Model
n = 20;
mode_choice = 1:3;
m = 1;
k = 1;
c = 1;
g = 0.5;

temp = spdiags(ones(n,1)*[-1 2 -1], [-1 0 1], sparse(n,n));
temp(1,1) = 1;

C = c * temp;
K = k * temp;
M = m*speye(n,n);

%% Nonlinearity
S = @(x)g*SP_nonlinearity(x);
DS = @(x)g*SP_nonlinearity_derv(x);

amplitudes = [1 2 4 8 ];    % amplitudes for which FRCs are computed
ncontsteps = 100;           % number of continuation steps

%% Steady-State tool
System.M=M;System.C=C;System.K=K;System.sys_order='second';
SSfull = SSR(System);               % full system
SSred = SSR(System,mode_choice);    % reduced system

SSred.S = S;
SSred.DS = DS;
SSfull.S = S;
SSfull.DS = DS;

%% Sequential continuation
omega1 = min(SSred.omega);
Omega_range = [0.7 1.3]*omega1;        
T_range = 2*pi./Omega_range;

for i = 1:length(amplitudes)
    alpha = amplitudes(i)*10^(-5);
    f0 = ones(n,1);                     % loading shape
    f = @(t,T)alpha*f0*sin(2*pi*t/T);   % loading function
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
legend('original full','reformulated full','reformulated reduced','location','SE')
xlabel('excitation frequency $$\Omega$$');
ylabel('response $$||x||_2$$')

xlim(Omega_range)
annotation('textbox', [0.13 0.275 .07 .05],'String',' $$F=0.01$$ ','LineStyle','none','FontSize',14,'interpreter','latex');
annotation('textbox', [0.13 0.41 .07 .05],'String',' $$F=0.02$$ ','LineStyle','none','FontSize',14,'interpreter','latex');
annotation('textbox', [0.13 0.54 .07 .05],'String',' $$F=0.04$$ ','LineStyle','none','FontSize',14,'interpreter','latex');
annotation('textbox', [0.13 0.67 .07 .05],'String',' $$F=0.08$$ ','LineStyle','none','FontSize',14,'interpreter','latex');
