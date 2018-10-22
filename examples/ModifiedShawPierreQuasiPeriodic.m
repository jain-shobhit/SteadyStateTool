%% Modified Shaw and Pierre example with quasiperiodic forcing 
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

SS = SSR(M,C,K);    % Instantiating the SSR package
SS.domain = 'freq';

%% External forcing
alpha = 0.01;   % Loading amplitude
f =  @(t,Omega) alpha*[1; 1]*sin(Omega*t);
f_theta = @(theta) alpha*[1; 1]*sin(2*pi*theta); 

%% single frequency linear
Omega = 1; % single frequency
SS.Omega = Omega;
SS.f_theta = f_theta;
SS.nh = 1; % number of harmonics to be used
SS.m = 50; % number of discretiztion along each coordinate on the Torus

[x_lin_freq, xd_lin_freq] = SS.LinearResponse();

% verification of the response computed using integral equations
T = 2*pi./Omega;
A = [zeros(n,n) eye(n,n);
    -M\K    -M\C];
Gamma = @(t,Omega)[zeros(n,1); M\f(t,Omega)]; % loading function
G = @(z)[ zeros(n,1); -M\S(z(1:n)) ];
F = @(t,z,Omega) A*z + G(z) + Gamma(t,Omega);
F_lin = @(t,z,Omega) A*z + Gamma(t,Omega);
% linear solution
[~, z_lin] = ode45(@(t,z)F_lin(t,z,Omega), [0 100*T], zeros(2*n,1)); % Transients
[t, z_lin] = ode45(@(t,z)F_lin(t,z,Omega), [0 3*T], z_lin(end,:)'); % Approximate periodic orbit

figure;
plot(t,z_lin(:,1:2),'DisplayName','ode45')
hold on;plot(2*pi*SS.theta_set/Omega,x_lin_freq, 'Displayname', 'Integral Eq.' ) 
grid on; axis tight

% %% multi-frequency, nonlinear
% S = @(x)[g * x(1)^3; 0];
% DS = @(x)[3 * g * x(1)^2, 0 ; 0, 0];
% SS.S = S;           % Setting the nonlinearity
% SS.Omega = [1, pi]; % Multifrequency
% SS.nh = 3;  % change number of along each base frequency to capture nonlinearity
% 
% % Nonlinear response with Picard iteration
% [x_picard,xd_picard] = SS.Picard('init',x_lin_freq,'tol',1e-3,'maxiter',20); % 'init', 'maxiter', 'tol' are optional arguments
% t_picard = SS.t;
% 
% % Nonlinear response with Newton--Raphson iteration
% SS.DS = DS;         % Setting the derivative
% [x_newton,xd_newton] = SS.NewtonRaphson(); % find out response using Picard iteration
% t_newton = SS.t;
