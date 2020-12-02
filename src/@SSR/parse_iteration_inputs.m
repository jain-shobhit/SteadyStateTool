function [x0,tol,maxiter] = parse_iteration_inputs(O,inputs)
% this method used for parsing interation inputs
defaultTol = 1e-6;
defaultInit = O.LinearResponse();
defaultmaxiter = 20;
p = inputParser;
addParameter(p,'tol',defaultTol)
addParameter(p,'init',defaultInit)
addParameter(p,'maxiter',defaultmaxiter)
parse(p,inputs{:});
x0 = p.Results.init;
maxiter = p.Results.maxiter;
tol = p.Results.tol;
end