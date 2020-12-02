function [pars] = parameters()
%% Parameters
pars.units = 'mm';  % whether to use SI units or mm units.
pars.L = 20; pars.H = 40;         % Dimensions
pars.n = 6;%6                      % No. of elements along length
pars.BC = 'B';                      % 'B' - simply supported on two opposite sides.
pars.t =  0.8;                       % 'DC' - clamped along one edge

pars.w= pars.H/20;    % curvature parameter (typical value = H/20)
if pars.w ~=0
    pars.R = 1/2*(pars.w^2+(pars.H/2)^2)/pars.w;
end

pars.p = 1; %20                       % Uniform Normal Pressure Amplitude
if strcmp(pars.BC, 'DC')
    pars.outcoord = [pars.L/2 pars.H/2 pars.w]; %[pars.L/2 0 pars.w];
else
    pars.outcoord = [pars.L/2 pars.H/2 pars.w];     % coordinates of output node
end
pars.lambda_dyn = 1 ;   % 0.01         % Load prefactor for time integration
pars.LinEq = 1;                         % Whether to solve for Linear Eqm
pars.modes = 1:5;
pars.M = length(pars.modes);                            % number of VMs
pars.PlotVM = [1 2];               % which VMs to plot
pars.Damping = 1;                       % 0 - No Damping, 1 - Rayleigh Damping
pars.dampc = 0.04;                      % damping coefficient for Rayleigh Damping
disp(pars)
