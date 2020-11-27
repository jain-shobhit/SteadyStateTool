function [Ksens]=Assemble_tangent_stiffness_derv(model,V)
%


ndof=model(1).ndof;
nelem=model(1).nelem;


M = size(V,2);          % No. of Modes in basis
Ksens = cell(1,M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for p = 1:M
    q = V(:,p); % give displacement in direction of pth mode
    
    I = zeros(nelem*18*18,1);
    J = zeros(nelem*18*18,1);
    KK = zeros(nelem*18*18,1);
    sh = 0;
    for i=1:nelem
        
        % extract element displacement
        ue=q(model(i).index);
        % rotate the displacement vector in local coordinates
        qe=model(i).T*ue;
        %     elnodes = model(i).node;
        
        [Ktel]=Allman_tangent_stiffness_derv(model(i).BL,model(i).Kxx,model(i).Kyy,model(i).Kxy,model(i).Kt,model(i).Am,model(i).area,qe);
        
        Ktel=model(i).T'*Ktel*model(i).T;
        
        
        for col = 1:18
            
            I(sh+18*(col-1)+1:sh+18*col,1) = ones(18,1)*model(i).index(col);
            J(sh+18*(col-1)+1:sh+18*col,1) = model(i).index;
            KK(sh+18*(col-1)+1:sh+18*col,1) = Ktel(col,:);
        end
        
        sh = sh + 324;
        
        
    end
    
    
    K = sparse(I,J,KK,ndof,ndof);
    
    Ksens{p} = K;
    clear q K
end