function S = Nonlinearity(model,Misc,x)
q = zeros(Misc.ndof,1);
q(Misc.freedofs) = x;
[~,F] = Assemble_Tangent_Stiffness(model,q);
S = F - Misc.Km * q;
S = S(Misc.freedofs);
end
