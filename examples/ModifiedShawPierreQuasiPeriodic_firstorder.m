%% Modified Shaw and Pierre example with quasiperiodic forcing 
% m*x1dd + k(x1-x2) + c(x1d-x2d) + g*x1^3 = f1(t)
% m*x2dd + k(2*x2-x1) + c(2*x2d - x1d) = f2(t)

clear all
%clc; close all

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

A = [ M , zeros(n) ; zeros(n), -K];
B = [ zeros(n) , M ; M       , C ];

System.A = A;
System.B = B;
System.sys_order = 'first';
SS = SSR(System);    % Instantiating the SSR package
SS.type = 'qp';     % Solve for quasiperiodic response (i.e. Fourier domain)

%% External forcing
alpha = 0.01;   % Loading amplitude
f =  @(t,Omega) alpha*[1; 1]*sin(Omega*t);
f_theta = @(theta) alpha*[1; 1]*sin(2*pi*theta); 
F_theta = @(theta) [ zeros(n,numel(theta));f_theta(theta)];

%% Linear, periodic forcing
Omega = 0.1; % single frequency
SS.Omega = Omega;
SS.F_theta = F_theta;
SS.nh = 1; % number of harmonics to be used
SS.m = 50; % number of discretiztion along each coordinate on the Torus


[z_lin_freq] = SS.LinearResponse();

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
hold on;plot(2*pi*SS.theta_set/Omega,z_lin_freq(n+1:end,:), 'Displayname', 'Integral Eq.' ) 
grid on; axis tight
legend('ode45', 'Integral Eq.')

%% Nonlinear, Quasi-periodic forcing 
S = @(x)[g * x(n+1)^3; 0];
DS = @(x)[3 * g * x(n+1)^2, 0 ; 0, 0];

R  = @(x) [zeros(n,1); S(x)];
DR = @(x) [zeros(n,2*n); zeros(n), DS(x)];

SS.R = R;           % Setting the nonlinearity
SS.DR = DR;         % Setting the derivative

f_theta = @(theta) alpha*[1 1; 0 0]*sin(2*pi*theta); 
F_theta = @(theta) [ zeros(n,size(theta,2));f_theta(theta)];
SS.F_theta = F_theta;
SS.nh = 3;          % change number of along each base frequency to capture 
                    % nonlinearity
                    
SS.Omega = [1,pi]; % Multifrequency                   
[z_lin_freq] = SS.LinearResponse();

% Nonlinear response with Picard iteration
try
[z_picard] = SS.Picard('init',z_lin_freq,'tol',1e-6,'maxiter',20); % 'init', 'maxiter', 'tol' are optional arguments
t_picard = SS.t;norm([z_picard],'fro')/size(SS.theta_set,2)
catch
end

% Nonlinear response with Newton--Raphson iteration
[z_newton] = SS.NewtonRaphson('maxiter',20); % find out response using Picard iteration
t_newton = SS.t;

%% sequential continuation to obtain response surface
param = 1;
n_cont_steps = 8;

Omega1_range=[0.9, 1.05];
Omega2=linspace(sqrt(3)-0.1,sqrt(3)+0.1,param);

Resp = nan(param,2*n_cont_steps+1);
for j = 1:length(Omega2)
    Omega0 = [Omega1_range(1), Omega2(j)];
    [Omega1_array,cont_sol,picard_flag] = SS.sequential_continuation(Omega1_range,'idx',1,'ncontsteps',n_cont_steps,'Omega0',Omega0);    
    % compute norm of the solution
    N = nan(length(Omega1_array),1);
    for k = 1:length(Omega1_array)
        if ~isempty(cont_sol{k})
            N(k) = norm(cont_sol{k},'fro')/size(SS.theta_set,2);
        end
    end
    Resp(j,:) = N;
end


figure; surf(Omega1_array, Omega2,Resp)
xlabel('$$\Omega_1$$ [rad/s] ')
ylabel('$$\Omega_2$$ [rad/s] ')
zlabel('$$\|\mathbf{z}\|_2$$ ')

