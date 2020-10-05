%% Oscillator chain with coupled nonlinearities
% This example is described in Section 4.2 (Figure 6) of the following 
% article 
% S. Jain, T. Breunung, G. Haller, Fast Computation of Steady-State 
% Response for Nonlinear Vibrations of High-Degree-of-Freedom Systems, 
% Nonlinear Dyn (2019) 97: 313. https://doi.org/10.1007/s11071-019-04971-1

clear all
clc


%%
n = 20;
m = 1;
k = 1;
c = 1;
g = 0.5;

temp = spdiags(ones(n,1)*[-1 2 -1], [-1 0 1], sparse(n,n));

C = c * temp;
K = k * temp;
M = m*speye(n,n);


S = @(x)g*SP_nonlinearity(x);
DS = @(x)g*SP_nonlinearity_derv(x);

alpha = 5e-5;
f0 = ones(n,1); % loading shape
f = @(t,T)alpha*f0*sin(2*pi*t/T); % loading function

%% Linear response
System.M=M;System.C=C;System.K=K;System.sys_order='second';
SS = SSR(System);    % Instantiating the SSR package
omega1 = SS.omega(1);
Omega_range = [0.7 1.3]*omega1;
T_range = 2*pi./Omega_range;
T = mean(T_range);
SS.f = f;           % setting the forcing
SS.T = T;           % setting the time period at which periodic loading is applied
SS.n_steps = 50;    % setting the number of time intervals within the time period
SS.order = 1;       % setting the order of integration, permissible values
                    % are 0,1,2,3,4

[x_lin, xd_lin] = SS.LinearResponse();
t_lin = SS.t;

%% Nonlinear response
SS.S = S;           % Setting the nonlinearity
% Nonlinear response with Picard iteration
[x_picard,xd_picard] = SS.Picard('init',x_lin,'tol',1e-6,'maxiter',20); % 'init', 'maxiter', 'tol' are optional arguments
t_picard = SS.t;
% hold on; plot(SS.t, x_picard, 'DisplayName', 'Picard');

% Nonlinear response with Newton--Raphson iteration
SS.DS = DS;         % Setting the derivative
[x_newton,xd_newton] = SS.NewtonRaphson(); % find out response using Picard iteration
t_newton = SS.t;
% hold on; plot(SS.t, x_newton, 'DisplayName', 'Newton');


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


