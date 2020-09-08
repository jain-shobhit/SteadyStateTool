function [S] = NonlinearityVK(model,x)
% The finite element model is taken from the following article: 
% Jain, S., Tiso, P., & Haller, G. (2018). Exact nonlinear model reduction 
% for a von Kármán beam: slow-fast decomposition and spectral submanifolds. 
% Journal of Sound and Vibration, 423, 195–211. 
% https://doi.org/10.1016/J.JSV.2018.01.049

q = zeros(model.nDOF,1);
q(model.freeDOFs) = x;
[~,F] = AssembleTangentStiffness(model,q);
S = F - model.K * q;
S = S(model.freeDOFs);
end

