function M = AssembleMass(model)
% The finite element model is taken from the following article: 
% Jain, S., Tiso, P., & Haller, G. (2018). Exact nonlinear model reduction 
% for a von Kármán beam: slow-fast decomposition and spectral submanifolds. 
% Journal of Sound and Vibration, 423, 195–211. 
% https://doi.org/10.1016/J.JSV.2018.01.049

ndof = model.nDOF;
nelem = model.nElements;
I = zeros(nelem*6*6,1);
J = zeros(nelem*6*6,1);
KK = zeros(nelem*6*6,1);
el = 0;
for i = 1:model.nElements
    index = 3*(i-1)+1:3*(i+1);
    
    Te = compute_rotation_matrix(model.nodes(i:i+1,:));
    
    Mel = ElementMass(model);
    Mel = Te.' * Mel * Te; 
        
    for col = 1:6
        I(el+6*(col-1)+1:el+6*col,1) = ones(6,1)*index(col);
        J(el+6*(col-1)+1:el+6*col,1) = index;
        KK(el+6*(col-1)+1:el+6*col,1) = Mel(col,:);
    end
    el = el + 36; 
  
end
M = sparse(I,J,KK,ndof,ndof);
