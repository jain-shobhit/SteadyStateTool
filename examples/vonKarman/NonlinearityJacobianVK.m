function [DS] = NonlinearityJacobianVK(model,x)

q = zeros(model.nDOF,1);
q(model.freeDOFs) = x;

[Kt,~] = AssembleTangentStiffness(model,q);

DS = Kt - model.K;
DS = DS(model.freeDOFs,model.freeDOFs);
end

