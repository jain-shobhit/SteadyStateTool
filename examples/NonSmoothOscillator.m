%% 2-DOF oscillator with non smooth nonlinearity
clear all
clc; close all

%% Properties
m = 1;
k = 1;
c = 0.02;
g = 0.5;
%% Linear Matrices
C = c*[2 -1;
    -1 2];
K = k* [2 -1;
    -1 2];
M = [m 0;
    0 m];
n = length(K);

%% External forcing
alpha = 0.02;   % Loading amplitude
f0 = [1;0];     % loading shape
f = @(t,T)alpha*f0*sin(2*pi*t/T); % loading function handle (also takes
% loading time period T as argument)

%% Linear response
Omega = 0.5;       % forcing frequency (rad/sec)
T = 2*pi/Omega;     % forcing time period
SS = SSR(M,C,K);    % Instantiating the SSR package
SS.f = f;           % setting the forcing
SS.T = T;           % setting the time period at which periodic loading is applied
SS.n_steps = 50;    % setting the number of time intervals within the time period
SS.order = 1;       % setting the order of integration, permissible values
                    % are 0,1,2,3,4

[x_lin, xd_lin] = SS.LinearResponse();
t_lin = SS.t;

%% nonsmooth oscillation
alp = 0.1; bet = 0.1;
S = @(x)[nonsmooth_nonlinearity(x(1),alp,bet); 0];
DS = @(x)[nonsmooth_nonlinearity_derv(x(1),alp,bet), 0 ; 0, 0];

SS.S = S;           
SS.DS = DS;

% Nonlinear response with Picard iteration
[x_picard,xd_picard] = SS.Picard(); % 'init', 'maxiter', 'tol' are optional arguments
t_picard = SS.t;
% hold on; plot(SS.t, x_picard, 'DisplayName', 'Picard');

% Nonlinear response with Newton--Raphson iteration
SS.DS = DS;         % Setting the derivative
[x_newton,xd_newton] = SS.NewtonRaphson(); % find out response using Picard iteration
t_newton = SS.t;
% hold on; plot(SS.t, x_newton, 'DisplayName', 'Newton');

%% sequential continuation
Omega_range = [0.5 1.5];
T_range = 2*pi./Omega_range;
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