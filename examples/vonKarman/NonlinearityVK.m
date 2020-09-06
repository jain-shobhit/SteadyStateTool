function [S] = NonlinearityVK(model,x)
q = zeros(model.nDOF,1);
q(model.freeDOFs) = x;
[~,F] = AssembleTangentStiffness(model,q);
S = F - model.K * q;
S = S(model.freeDOFs);
end

