function [K,F] = AssembleTangentStiffness(model,u)

ndof = model.nDOF;
nelem = model.nElements;
I = zeros(nelem*6*6,1);
J = zeros(nelem*6*6,1);
KK = zeros(nelem*6*6,1);

el = 0;
F = zeros(ndof,1);
for i = 1:model.nElements
    index = 3*(i-1)+1:3*(i+1);
    
    Te = compute_rotation_matrix(model.nodes(i:i+1,:));
    
    qe = Te * u(index); 
    [Kel,Fel] = ElementStiffness(model,qe);
    
    Kel = Te.' * Kel * Te; 
    F(index) = F(index)+ Te.' * Fel;
    
    for col = 1:6
        I(el+6*(col-1)+1:el+6*col,1) = ones(6,1)*index(col);
        J(el+6*(col-1)+1:el+6*col,1) = index;
        KK(el+6*(col-1)+1:el+6*col,1) = Kel(col,:);
    end
    el = el + 36; 
  
end
K = sparse(I,J,KK,ndof,ndof);
