function [x0,eta0,tol,maxiter] = parse_iteration_inputs_ref(O,inputs)
% this method used for parsing interation inputs
p = inputParser;
defaultTol = 1e-6;
a = 1;
[x0,~,defaultEta] = O.LinearResponse();

% check if user specified init
for j = 1:length(inputs)
    tf = strcmp(inputs{j},'init');
    if tf 
        a = 0;
        break
    end
end

% if user specified init, avoid the unnecessary function evaluation
if a    
    defaultInit = - O.U.' * SSR.evaluate_fun_over_array(O.S,x0,false);
    addParameter(p,'init',defaultInit)    
else
    addParameter(p,'init',0)
end

defaultmaxiter = 20;
addParameter(p,'eta0',defaultEta)
addParameter(p,'tol',defaultTol)
addParameter(p,'maxiter',defaultmaxiter)
parse(p,inputs{:});
x0 = p.Results.init;
maxiter = p.Results.maxiter;
tol = p.Results.tol;
eta0 = p.Results.eta0;
end