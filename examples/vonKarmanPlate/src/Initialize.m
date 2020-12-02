function [model,Misc,pars] = Initialize(pars)
% model generation

[Misc, model] = PlateModel(pars);

% stiffness matrix assembly
disp('Assembling the material stiffness matrix at zero displacement...')
% tic
[Misc.Km,model]=Assemble_Linear_Stiffness(model);
Misc.KmAssemblytime=toc;
% mass matrix assembly
disp('Assembling the mass matrix...')
% tic
[Misc.Mm,model]=Assembly_M(Misc.nodes,Misc.elements,model);
Misc.MmAssemblytime=toc;


%% preload perfect structure
if pars.LinEq
disp('Solving for the prebuckling solution...')
% tic
Misc.qL = zeros(Misc.ndof,1);
Misc.qL(Misc.freedofs,1)  = ...
    Misc.Km(Misc.freedofs,Misc.freedofs)\Misc.loads(Misc.freedofs);
% toc
outdof = [Misc.outnode*6-4 Misc.outnode*6-5 Misc.outnode*6-3];
u_out = Misc.qL(outdof);
disp(['u_out_lin = ' num2str(u_out(2,end)) ', v_out_lin = ' num2str(u_out(1,end)) ', w_out_lin = ' num2str(u_out(3,end))])

% figure
% def_plot(Misc.qL,0.1,Misc.nodes,Misc.elements)
% title('Linear Equilibrium Solution')
end

%% Free Vibration
%  vibration modes @ load = 0
disp('Solving the vibration problem at load = 0 ...')

% tic
Misc.Mmred = Misc.Mm(Misc.freedofs,Misc.freedofs);
Misc.Kmred = Misc.Km(Misc.freedofs,Misc.freedofs);
[Misc.V,Misc.D] = eigs(Misc.Kmred,Misc.Mmred,max(pars.modes),'SM');   % Return M smallest eigenvalues and vectors
[Misc.omega2,Misc.V] = sortVMs(Misc.V,Misc.D);
Misc.VMs = zeros(Misc.ndof,pars.M);
Misc.VMs(Misc.freedofs,:) = Misc.V(:,pars.modes);
Misc.EigenValueProblemSolutiontTime = toc;

% unit mass scaing of VMs
Misc.mu = diag(Misc.VMs'*Misc.Mm*Misc.VMs);
Misc.VMs(Misc.freedofs,:) = Misc.VMs(Misc.freedofs,:)*diag(1./sqrt(Misc.mu));

% if ~isempty(pars.PlotVM)
%     for i = 1:length(pars.PlotVM)
%         figure
%         def_plot(Misc.VMs(:,pars.PlotVM(i)),0.1,Misc.nodes,Misc.elements)
%         title(strcat('Vibration Mode ',num2str(pars.PlotVM(i)),...
%             '  frequency = ',num2str(sqrt(Misc.omega2(i))/(2*pi)),'  Hz') )
%     end
% end
%% Damping
% Rayleigh Damping
if pars.Damping == 1 
    disp('Rayleigh damping Selected')
    W =   sqrt(Misc.omega2(1:2));
    Misc.a = [W(1) 1/W(1);W(2) 1/W(2)]\[pars.dampc; pars.dampc];
    Misc.Cm = Misc.a(2)*Misc.Mm + Misc.a(1)*Misc.Km;    
    for i = 1:model(1).nelem
        model(i).Cel = Misc.a(2)*model(i).Mel + Misc.a(1)*model(i).Kel;
    end
else    
    disp('Undamped system')
    Misc.Cm = sparse(model(1).ndof,model(1).ndof);    
end

end

