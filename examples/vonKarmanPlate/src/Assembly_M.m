function [M,model]=Assembly_M(nodes,elements,model)

ndof=model(1).ndof;
nelem=size(elements,1);
I = zeros(nelem*18*18,1);
J = zeros(nelem*18*18,1);
KK = zeros(nelem*18*18,1);
sh = 0;

tic
for i=1:nelem
    
    node1=[nodes(elements(i,2),2) nodes(elements(i,2),3) nodes(elements(i,2),4)];
    node2=[nodes(elements(i,3),2) nodes(elements(i,3),3) nodes(elements(i,3),4)];
    node3=[nodes(elements(i,4),2) nodes(elements(i,4),3) nodes(elements(i,4),4)];
    
    %     x=[node1(1) node2(1) node3(1)];
    %     y=[node1(2) node2(2) node3(2)];
    
    elnodes=[node1; node2; node3];
    
    T=rot3D(elnodes);
    %rotate element coordinates into local system
    elnodes=elnodes*T(1:3,1:3)';
    
    
    t = model(i).property(4);
    rho = model(i).property(5);
    
    % membrane contribution
    [Melm,~]=Allman_membrane_mass(elnodes,t,rho);
    
    Melmex=zeros(18,18);
    Melmex([1 2 6 7 8 12 13 14 18],[1 2 6 7 8 12 13 14 18])=Melm;
    Melmex=T'*Melmex*T;
    
    % bending contribution
    [Melb,~]=Allman_bending_mass(elnodes,t,rho);
    
    Melbex=zeros(18,18);
    Melbex([3 4 5 9 10 11 15 16 17],[3 4 5 9 10 11 15 16 17])=Melb;
    Melbex=T'*Melbex*T;
    Mel=Melbex+Melmex;
    model(i).Mel = Mel;
    
    % % consistent mass matrix
    %      M(index,index)=M(index,index)+Mel;
    
    for col = 1:18
        I(sh+18*(col-1)+1:sh+18*col,1) = ones(18,1)*model(i).index(col);
        J(sh+18*(col-1)+1:sh+18*col,1) = model(i).index;
        KK(sh+18*(col-1)+1:sh+18*col,1) = Mel(col,:);
    end
    sh = sh + 324;
    
end

M = sparse(I,J,KK,ndof,ndof);

