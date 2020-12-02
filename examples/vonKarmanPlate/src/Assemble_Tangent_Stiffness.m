function [K,Fint]=Assemble_Tangent_Stiffness(model,u)
%
% [K,Fint]=Assemble_shell_new(model,q)
%
% K = tangential stiffness matrix
% Fint = internal forces

ndof=model(1).ndof;
nelem=model(1).nelem;
Fint=zeros(ndof,1);
I = zeros(nelem*18*18,1);
J = zeros(nelem*18*18,1);
KK = zeros(nelem*18*18,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sh = 0;

for i=1:nelem
    
    % extract element displacement
    ue=u(model(i).index);
    % rotate the displacement vector in local coordinates
    qe=model(i).T*ue;
    % Calculate element tangent stiffness and Internal forces
    [Fel,Ktel]=Allman_tangent_stiffness_new(model(i).BL, ...
        model(i).Kxx,model(i).Kyy,model(i).Kxy,model(i).Kt,...
        model(i).Am,model(i).area,qe);
    % Rotate K and F to global coordinates
    Ktel=model(i).T'*Ktel*model(i).T;    
    Fel=model(i).T'*Fel;
    
    % scatter and assemble    
    for col = 1:18
        
        I(sh+18*(col-1)+1:sh+18*col,1) = ones(18,1)*model(i).index(col);
        J(sh+18*(col-1)+1:sh+18*col,1) = model(i).index;
        KK(sh+18*(col-1)+1:sh+18*col,1) = Ktel(col,:);
    end    
    sh = sh + 324;  
    
    Fint(model(i).index)=Fint(model(i).index)+Fel;
end



K = sparse(I,J,KK,ndof,ndof);