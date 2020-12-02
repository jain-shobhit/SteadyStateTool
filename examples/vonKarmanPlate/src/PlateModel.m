function [Misc, model] = PlateModel(pars)
% This function creates the model for plate using the given parameters in
% the pars data structure
pars.m = ceil(pars.H/pars.L*pars.n) ;
[Misc.nodes,Misc.elements,~,Misc.freedofs,Misc.prop]=cantilever_2(pars.L,pars.H,pars.n,...
    pars.m,pars.BC,pars.w,pars.t);

Misc.ndof=size(Misc.nodes,1)*6;
% creating the FE database
model = model_database_2(Misc.nodes,Misc.elements,Misc.prop);

% uniform extrenal pressure
[Misc.loads]=Assembly_pressure(model,pars.p);
model(1).nbc_el = 1:model(1).nelem;

Misc.alldofs = 1:Misc.ndof;
Misc.bounddofs = setdiff(Misc.alldofs,Misc.freedofs);
model(1).nbc_el = 1:model(1).nelem;
dist = zeros(size(Misc.nodes,1),1);
for j = 1:size(Misc.nodes,1)
    dist(j) = norm(Misc.nodes(j,2:4) - pars.outcoord); 
end

[~,outnode] = min(dist);
Misc.outnode = outnode;
Misc.out_dof = Misc.outnode*6-3; % z direction displacement
