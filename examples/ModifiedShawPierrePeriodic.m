%% Modified Shaw and Pierre example: replace spring 1 with nonlinear spring x1^3
% m*x1dd + k(x1-x2) + c(x1d-x2d) + g*x1^3 = f1(t)
% m*x2dd + k(2*x2-x1) + c(2*x2d - x1d) = f2(t)

% This example is described in Section 4.1 (Figures 1,2) of the following 
% article 
% S. Jain, T. Breunung, G. Haller, Fast Computation of Steady-State 
% Response for Nonlinear Vibrations of High-Degree-of-Freedom Systems, 
% Nonlinear Dyn (2019) 97: 313. https://doi.org/10.1007/s11071-019-04971-1


clear all
clc; close all

%% Properties
m = 1;
k = 1;
c = 0.3;
g = 0.5;
%% Linear Matrices
C = c*[1 -1;
    -1 2];
K = k* [1 -1;
    -1 2];
M = [m 0;
    0 m];
n = length(K);

%% External forcing
alpha = 2e-2;   % Loading amplitude
f0 = [1;1];     % loading shape
f = @(t,T)alpha*f0*sin(2*pi*t/T); % loading function handle (also takes
% loading time period T as argument)

%% Linear response
Omega = 0.25;       % forcing frequency (rad/sec)
T = 2*pi/Omega;     % forcing time period

System.M=M;System.C=C;System.K=K;System.sys_order='second';
SS = SSR(System);    % Instantiating the SSR package
SS.f = f;           % setting the forcing
SS.T = T;           % setting the time period at which periodic loading is applied
SS.n_steps = 100;    % setting the number of time intervals within the time period
SS.order = 1;       % setting the order of integration, permissible values
% are 0,1,2,3,4

[x_lin, xd_lin] = SS.LinearResponse();
t_lin = SS.t;

%% Nonlinearity function handle
S = @(x)[g * x(1)^3; 0];
DS = @(x)[3 * g * x(1)^2, 0 ; 0, 0];

%% Nonlinear response
SS.S = S;           % Setting the nonlinearity
% Nonlinear response with Picard iteration
[x_picard,xd_picard] = SS.Picard('init',x_lin,'tol',1e-3,'maxiter',20); % 'init', 'maxiter', 'tol' are optional arguments
t_picard = SS.t;
% hold on; plot(SS.t, x_picard, 'DisplayName', 'Picard');

% Nonlinear response with Newton--Raphson iteration
SS.DS = DS;         % Setting the derivative
[x_newton,xd_newton] = SS.NewtonRaphson(); % find out response using Picard iteration
t_newton = SS.t;
% hold on; plot(SS.t, x_newton, 'DisplayName', 'Newton');

%% verification of the response computed using integral equations
A = [zeros(n,n) eye(n,n);
    -M\K    -M\C];

Gamma = @(t,T)[zeros(n,1); M\f(t,T)]; % loading function
G = @(z)[ zeros(n,1); -M\S(z(1:n)) ];
F = @(t,z,T) A*z + G(z) + Gamma(t,T);
F_lin = @(t,z,T) A*z + Gamma(t,T);

figure
% linear solution
[~, z_lin] = ode45(@(t,z)F_lin(t,z,T), [0 100*T], [0; 0;0;0]); % Transients
[t, z_lin] = ode45(@(t,z)F_lin(t,z,T), [0 3*T], z_lin(end,:)'); % Approximate periodic orbit
for j = 1:n
    subplot(n,1,j)
    plot(t, z_lin(:,j), 'DisplayName', 'Linear - ode45');
    hold on; plot(t_lin, x_lin(j,:), 'DisplayName', 'Linear - IE');
    xlabel('$$t$$'); ylabel(['$$x_' num2str(j) '$$'])
    axis tight; grid on;legend('show')
end

[~, z] = ode45(@(t,z)F(t,z,T), [0 100*T], zeros(2*n,1)); % Transients
[t, z] = ode45(@(t,z)F(t,z,T), [0 3*T], z(end,:)'); % Approximate periodic orbit

% plots
figure;
for j = 1:n
    subplot(n,1,j)
    plot(t, z(:,j), 'DisplayName', 'ode45')
    hold on; plot(t_picard, x_picard(j,:), 'DisplayName', 'Picard');
    plot(T+t_newton, x_newton(j,:), 'DisplayName', 'Newton');
    xlabel('$$t$$'); ylabel(['$$x_' num2str(j) '$$'])
    axis tight; grid on; legend('show')
    
end

%% sequential continuation
T_range = [6 21];
[Omega_array,cont_sol,~] = SS.sequential_continuation(T_range,'ncontsteps',100);
% compute norm of the solution
N = nan(length(Omega_array),1);
for j = 1:length(Omega_array)
    if ~isempty(cont_sol{j})
        N(j) = SS.dt*norm(cont_sol{j},'fro');
    end
end
figure; semilogy(Omega_array,N,'-b');axis tight; grid on;
xlabel('$$\Omega$$ [rad/s]'); ylabel('$$||\mathbf{x}||_2$$')

%% Ampltude variation
alpha_array = [0.01 0.02 0.04 0.08]; % Amplitude array
figure;
for k = 1:length(alpha_array)
    f = @(t,T)alpha_array(k)*f0*sin(2*pi*t/T); % loading function handle
    SS.f = f;
    [Omega_array,cont_sol,picard_flag] = SS.sequential_continuation(T_range,'ncontsteps',100);
    N = nan(length(Omega_array),1);
    for j = 1:length(Omega_array)
        if ~isempty(cont_sol{j})
            N(j) = SS.dt*norm(cont_sol{j},'fro');
        end
    end
    % plot Picard and Newton solutions in different color
    picard_idx = find(picard_flag);
    newton_idx = find(~picard_flag);
    N_picard = N;
    N_picard(newton_idx)=nan;    
    semilogy(Omega_array,N_picard,'-b','DisplayName','Picard');
    text(Omega_array(1),0.9*N(1),['$$A=$$' num2str(alpha_array(k))])
    hold on; 
    if ~isempty(newton_idx)
        newton_idx = [newton_idx(1)-1; newton_idx; newton_idx(end)+1];
        semilogy(Omega_array(newton_idx),N(newton_idx),'-r','DisplayName', 'Newton-Raphson');
    end
end
axis tight; grid on;
xlabel('$$\Omega$$ [rad/s]'); ylabel('$$||\mathbf{x}||_2$$')




