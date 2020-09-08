function [Misc,model] = Beam_Model(Geo,nElements,B)
% This finite element model is taken from the following article: 
% Jain, S., Tiso, P., & Haller, G. (2018). Exact nonlinear model reduction 
% for a von Kármán beam: slow-fast decomposition and spectral submanifolds. 
% Journal of Sound and Vibration, 423, 195–211. 
% https://doi.org/10.1016/J.JSV.2018.01.049

Damping = 2;    % Structural damping : 0 - Nil,  1 - Rayleigh, 2 - Linear viscoelastic

% Assumed a curved beam with N equidistant elements along the x axis
model.nElements = nElements;
model.nNodes = nElements +1;
model.nDOF = model.nNodes*3;
model.dx = Geo.l/model.nElements;
model.b = Geo.b;

model.nodes = [(0:model.dx:Geo.l).' zeros(model.nNodes,1)];
if Geo.a~=0
    R = 1/2*(Geo.a^2+(Geo.l/2)^2)/Geo.a;
    theta0 = asin(Geo.l/2/R);
    
    for i = 1:size(model.nodes,1)
        th = asin((model.nodes(i,1)-Geo.l/2)/R);
        model.nodes(i,2) = R *  cos(th) - R * cos(theta0);
    end
end

% Material Properties (Aluminium)
model.E = 70e9;     % Youngs Modulus (Pa)
model.A = Geo.A;    % Cross section area of Beam
model.I = Geo.I;    % Moment of inertia of beam
model.rho = 2700;   % Density
model.kappa = 1e8;  % Material damping for Kelvin Voigt model
     

% Loading
model.pw = 1;      % loading factor along the transverse direction
model.pu = 0;      % loading factor along the axial direction

model.loads = sparse(model.nDOF,1);
l = ones(model.nNodes,1);
l(2:end-1) = 2;
model.uDOF = 1:3:model.nDOF; % axial degrees of freedom
model.wDOF = setdiff(1:model.nDOF,model.uDOF); % transverse DOFs 
model.loads(2:3:end) = model.pw*(model.dx/2)*l; % 
model.loads(model.uDOF) = model.pu*(model.dx/2)*l;

% Boundary conditions
model.fixDOFs = BoundaryCoundition(model,B);
model.freeDOFs = setdiff(1:model.nDOF, model.fixDOFs);
model.uDOF = setdiff(model.uDOF,model.fixDOFs);
model.wDOF = setdiff(model.wDOF,model.fixDOFs);

% Stiffness
model.K = AssembleLinearStiffness(model);
model.M = AssembleMass(model);

%
model.modes = 1:5;
model.m = length(model.modes);
% Vibration problem
Mred = model.M(model.freeDOFs,model.freeDOFs);
model.Ks = model.K;
Kred = model.K(model.freeDOFs,model.freeDOFs);
[model.V,D] = eigs(Kred,Mred,max(model.modes),'SM');   % Return M smallest eigenvalues and vectors
[omega2,V] = sortVMsVK(model.V,Mred,D);
disp('First 5 undamped Natural Frequencies [Hz]')
disp(sqrt(omega2)/(2*pi))
model.omega2 = omega2(model.modes);
model.VMs = zeros(model.nDOF,model.m);
model.VMs(model.freeDOFs,:) = V(:,model.modes);
Misc.omega = sqrt(model.omega2);


figure
PlotDeformation(model,model.VMs(:,1),0.1);

switch Damping
    case 1
        % Rayleigh Damping
        W =   sqrt(model.omega2(1:2));
        a = [W(1) 1/W(1);W(2) 1/W(2)]\[0.02;0.02];
        model.C = (a(2)*model.M + a(1)*model.K);
               
    case 2
        % ViscoElastic Damping with linear strain rate
        model.C = (model.kappa/model.E)*model.Ks;        
end

