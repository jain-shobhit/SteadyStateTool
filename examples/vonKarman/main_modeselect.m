clear 
clc
close all

%% Von Karman Beam example
epsilon = 7e-3; % h/L 7e-3 1e-3
midp_height = 0.005; %0.005 height of midpoint relative to the ends (measure of curvature of the beam)
[Geo] = Geometry(epsilon,midp_height);
nElements = 10;
BC = 'B'; % simply supported
[Misc,model] = Beam_Model(Geo,nElements,BC);

n = length(model.freeDOFs);

%% Mass Matrix
M = model.M(model.freeDOFs,model.freeDOFs);

%% Stiffness Matrix
K = model.K(model.freeDOFs,model.freeDOFs);

%% Damping Matrix
C = model.C(model.freeDOFs,model.freeDOFs);

%% Transverse loading on Beam
A = 100;
Fext = @(t,T) A*model.loads(model.freeDOFs) * sin(2*pi*t/T); % specify forcing here if available

%% Steady-State tool
FULL = SSR(M,C,K,1:n);    
FULL.f = Fext;           
FULL.T = 2*pi/25;           
FULL.n_steps = 50;    
FULL.order = 1;       
                    
%% Mode selection 

% Input parameters for mode selection
param.curvtolerance = 0.05;
param.nmodes = 10;
param.type = 1; % set to 0 if param.nmodes is the max number of modes,
% set to 1 if param.nmodes is the desired number of modes
param.example = 'simple'; % always use simple if enough memory is available to store X below

X = S11(model); % this conains the quadratic nonlinearities for all dofs (X_i).
% X_i is symmetric, and can be found from the form M\ddot{q} + C\dot{q} +
% Kq +\Gamma(q) = F(t) (i denotes the specific row of this equation)
modes = modeselect(FULL,param,X)





