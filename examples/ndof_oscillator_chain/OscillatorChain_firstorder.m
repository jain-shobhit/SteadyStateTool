%% Oscillator chain with coupled nonlinearities
clear all
clc


%%
n = 20;
n_steps = 50;
m = 1;
k = 1;
c = 1;
g = 0.5;

temp = spdiags(ones(n,1)*[-1 2 -1], [-1 0 1], sparse(n,n));

C = c * temp;
K = k * temp;
M = m*speye(n,n);

A = [ M , zeros(n) ; zeros(n), -K];
B = [ zeros(n) , M ; M       , C ];

S = @(x)g*SP_nonlinearity_O1(x);
DS = @(x)g*SP_nonlinearity_derv_O1(x);

R  = @(x) [zeros(n,1); S(x)];
DR = @(x) [zeros(n,2*n); zeros(n), DS(x)];

alpha = 5e-5;
f0 = ones(n,1); % loading shape
f = @(t,T)alpha*f0*sin(2*pi*t/T); % loading function

F = @(t,T) [zeros(n,n_steps+1); f(t,T)];
%% Linear response
System.A = A;
System.B = B;
System.sys_order = 'first';
SS = SSR(System);    % Instantiating the SSR package

omega1 = imag(SS.lambda(1));
Omega_range = [0.7 1.3]*omega1;
T_range = 2*pi./Omega_range;
T = mean(T_range);
SS.F = F;           % setting the forcing
SS.T = T;           % setting the time period at which periodic loading is applied
SS.n_steps = n_steps;    % setting the number of time intervals within the time period
SS.order = 1;       % setting the order of integration, permissible values
                    % are 0,1,2,3,4

[z_lin] = SS.LinearResponse();
t_lin = SS.t;
%% Nonlinear response
SS.R = R;           % Setting the nonlinearity
% Nonlinear response with Picard iteration
[z_picard] = SS.Picard('init',z_lin,'tol',1e-6,'maxiter',20); % 'init', 'maxiter', 'tol' are optional arguments
t_picard = SS.t;
hold on; plot(SS.t, z_picard(n+1:end,:), 'DisplayName', 'Picard');

% Nonlinear response with Newton--Raphson iteration
SS.DR = DR;         % Setting the derivative
[z_newton] = SS.NewtonRaphson(); % find out response using Picard iteration
t_newton = SS.t;
%hold on; plot(SS.t, z_newton, 'DisplayName', 'Newton');

%% sequential continuation
tic
[Omega_array,cont_sol,~] = SS.sequential_continuation(T_range,'ncontsteps',100);
toc
% compute norm of the solution
N = nan(length(Omega_array),1);
for j = 1:length(Omega_array)
    if ~isempty(cont_sol{j})
        N(j) = SS.dt*norm(cont_sol{j},'fro');
    end
end
figure; semilogy(Omega_array,N,'-b');axis tight; grid on;
xlabel('$$\Omega$$ [rad/s]'); ylabel('$$||\mathbf{x}||_2$$')


