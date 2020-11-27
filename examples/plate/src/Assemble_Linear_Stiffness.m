function [K,model]=Assemble_Linear_Stiffness(model)

ndof=model(1).ndof;
nelem=model(1).nelem;
I = zeros(nelem*18*18,1);
J = zeros(nelem*18*18,1);
KK = zeros(nelem*18*18,1);
sh = 0;
tic
for i=1:nelem
    %rotate element coordinates into local system
    Kel = model(i).T'*model(i).Kt*model(i).T;  
    model(i).Kel = Kel;
    % scatter and Assemble
    for col = 1:18        
        I(sh+18*(col-1)+1:sh+18*col,1) = ones(18,1)*model(i).index(col);
        J(sh+18*(col-1)+1:sh+18*col,1) = model(i).index;
        KK(sh+18*(col-1)+1:sh+18*col,1) = Kel(col,:);
    end    
    sh = sh + 324;  

end

K = sparse(I,J,KK,ndof,ndof);
