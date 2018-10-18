%% Modified Shaw and Pierre example: replace spring 1 with nonlinear spring x1^3
% m*x1dd + k(x1-x2) + c(x1d-x2d) + g*x1^3 = f1(t)
% m*x2dd + k(2*x2-x1) + c(2*x2d - x1d) = f2(t)

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
alpha = 1e-1;   % Loading amplitude
f0 = [1;1];     % loading shape
f = @(t,T)alpha*f0*sin(2*pi*t/T); % loading function handle (also takes
% loading time period T as argument)
Df = @(t,T)-(alpha*2*pi/(T^2))*f0*(cos(2*pi*t/T).*t); % derivative of
% loading function
% w.r.t. T

%% Linear response
Omega = 0.25;       % forcing frequency (rad/sec)
T = 2*pi/Omega;     % forcing time period
SS = SSR(M,C,K);    % Instantiating the SSR package
SS.f = f;           % setting the forcing
SS.T = T;           % setting the time period at which periodic loading is applied
SS.n_steps = 50;    % setting the number of time intervals within the time period
SS.order = 1;       % setting the order of integration, permissible values 
                    % are 0,1,2,3,4

[x_lin, ~] = SS.LinearResponse();
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

[~, z] = ode45(@(t,z)F(t,z,T), [0 100*T], [0; 0;0;0]); % Transients
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
% T_range = [6 14];
% [T_array,cont_sol] = SS.sequential_continuation(T_range);


%% Quasi-periodic case
SS.domain = 'freq';
f_theta = @(theta,Omega) alpha*f0*sin(2*pi*theta); 
SS.Omega = Omega;
SS.f_theta = f_theta;
SS.nh = 10;
SS.m = 50;

[x_lin_freq, xd_lin_freq] = SS.LinearResponse();

figure;
plot(2*pi*SS.theta_set/Omega,x_lin_freq) 
hold on; plot(SS.t,x_lin); axis tight, grid on; 

[x_picard_freq,xd_picard_freq] = SS.Picard('init',x_lin_freq,'tol',1e-3,'maxiter',20); % 'init', 'maxiter', 'tol' are optional arguments

figure;
plot(2*pi*SS.theta_set/Omega,x_picard_freq) 
hold on; plot(t_picard,x_picard); axis tight, grid on; 


