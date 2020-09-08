function [DS] = NonlinearityJacobianVK(model,x)
% The finite element model is taken from the following article: 
% Jain, S., Tiso, P., & Haller, G. (2018). Exact nonlinear model reduction 
% for a von Kármán beam: slow-fast decomposition and spectral submanifolds. 
% Journal of Sound and Vibration, 423, 195–211. 
% https://doi.org/10.1016/J.JSV.2018.01.049

q = zeros(model.nDOF,1);
q(model.freeDOFs) = x;

[Kt,~] = AssembleTangentStiffness(model,q);

DS = Kt - model.K;
DS = DS(model.freeDOFs,model.freeDOFs);
end

