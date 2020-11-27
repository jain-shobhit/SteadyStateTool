function DS = NonlinearityJacobian(model,Misc,x)

q = zeros(Misc.ndof,1);
q(Misc.freedofs) = x;

[Kt,~] = Assemble_Tangent_Stiffness(model,q);

DS = Kt - Misc.Km;
DS = DS(Misc.freedofs,Misc.freedofs);
end
