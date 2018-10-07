%% Modified Shaw and Pierre example: replace spring 1 with nonlinear spring x1^3
% m*x1dd + k(x1-x2) + c(x1d-x2d) + g*x1^3 = f1(t)
% m*x2dd + k(2*x2-x1) + c(2*x2d - x1d) = f2(t)

clear all
clc

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
%% Nonlinearity function handle
S = @(x)[g * x(1)^3; 0];
DS = @(x)[3 * g * x(1)^2, 0 ; 0, 0];

%% External forcing 
alpha = 4e-2;   % Loading amplitude
f0 = [1;1];     % loading shape
f = @(t,T)alpha*f0*sin(2*pi*t/T); % loading function handle (also takes 
                                  % loading time period T as argument)
Df = @(t,T)-(alpha*2*pi/(T^2))*f0*(cos(2*pi*t/T).*t); % derivative of 
                                                      % loading function 
                                                      % w.r.t. T

%% Linear response
Omega = 0.25;   % forcing frequency

SS = SSR(M,C,K);    % Instantiating the SSR package    
SS.f = f;           % setting the forcing
SS.T = 2*pi/Omega;  % setting the time period at which periodic loading is applied 
SS.n_steps = 50;    % setting the number of time intervals within the time period

[x, xd] = SS.LinearResponse();
plot(SS.t, x); axis tight; grid on;

%% Nonlinear response
SS.S = S;           % Setting the nonlinearity
SS.DS = DS;         % Setting the derivative
