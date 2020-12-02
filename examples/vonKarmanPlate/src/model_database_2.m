function model = model_database_2(nodes,elements,prop)
% INPUT:
% nodes = nodal matrix [node_id x y z]
% elements = element connectivity matrix [el_id node1 node2 node3]
% prop = [prop_id E nu t rho]
%
% OUTPUT:
% model(i).T = rotation matrix of the element
% model(i).node = element nodal coordinates
% model(i).index = element connectivity
% model(i).Am = isotropic constitutive matrix for the membrane part
% model(i).area = Area of the element 
% model(i).xij = xij --> xi-xj
% model(i).yij = yij --> yi-yj 
% model(i).BL = linear displacement-strain matrix;
% model(i).Kxx = Kxx, quadratic displacement-strain matrix, (w,x)^2 = q'*Kxx*q;
% model(i).Kyy = Kyy, quadratic displacement-strain matrix, (w,y)^2 = q'*Kyy*q;
% model(i).Kxy = Kxy, quadratic displacement-strain matrix, w,x+w,y = q'*Kxy*q;
% N.B. eps = BL*q+1/2*BNL*q, where:
%      eps = strain vector, [ex ey exy kx ky kxy]'
%      q   = displacement vector, [ux uy uz rx ry rz]'
%      BNL = [q'*Kxx q'*Kyy q'*Kxy 0 0 0];
% model(i).Kt = shell stiffness matrix;

nelem = size(elements,1);
ndof = size(nodes,1)*6;

model.ndof = ndof;
model.nelem = nelem;

% for Allman's '88 membrane element (see Ref. page 14)
ab=1;
b0=4/9;
b = [1/12 5/12 1/2 0 1/3 -1/3 -1/12 -1/2 -5/12];

for i = 1:nelem
    
    node1=[nodes(elements(i,2),2) nodes(elements(i,2),3) nodes(elements(i,2),4)]; 
    node2=[nodes(elements(i,3),2) nodes(elements(i,3),3) nodes(elements(i,3),4)];
    node3=[nodes(elements(i,4),2) nodes(elements(i,4),3) nodes(elements(i,4),4)];
    
%     x=[node1(1) node2(1) node3(1)];
%     y=[node1(2) node2(2) node3(2)];
    
    elnodes=[node1; node2; node3];
    
    % element rotation matrix
    model(i).T = rot3D(elnodes);
    
    %element nodal coordinates
    model(i).node = elnodes*model(i).T(1:3,1:3)';
    
    conn=elements(i,:);
    index=[6*(conn(2)-1)+1:6*(conn(2)-1)+6 6*(conn(3)-1)+1:6*(conn(3)-1)+6 6*(conn(4)-1)+1:6*(conn(4)-1)+6  ];
    
    % element connectivity
    model(i).index = index;
    
    %element property
    model(i).property = prop(elements(i,5),:);
    
    E = prop(elements(i,5),2); nu = prop(elements(i,5),3); h = prop(elements(i,5),4);
    Am=h*E/(1-nu^2)*[ 1 nu        0;
                      nu 1        0;
                      0  0 (1-nu)/2];
                  
    model(i).Am = Am;
   
     % geometric quantities
     
    x=model(i).node(:,1); y=model(i).node(:,2);
     
    x12=x(1)-x(2); x23=x(2)-x(3); x31=x(3)-x(1); x21=-x12; x32=-x23; x13=-x31;
    y12=y(1)-y(2); y23=y(2)-y(3); y31=y(3)-y(1); y21=-y12; y32=-y23; y13=-y31;
    
    xij = [x12 x23 x31 x21 x32 x13];
    yij = [y12 y23 y31 y21 y32 y13];
    
    % element area
    A=(y21*x13-x21*y13)/2; A2=2*A; A4=4*A;
    
    model(i).area = A; 
    model(i).xij = xij;
    model(i).yij = yij;
    
    BL=sparse(3,18);
    BL(:,[1 2 6 7 8 12 13 14 18])= [[y23,0,x32];[0,x32,y23];[y23*(y13-y21),x32*(x31-x12),(x31*y13-x12*y21)*2]*ab/6; ...
                                [y31,0,x13];[0,x13,y31];[y31*(y21-y32),x13*(x12-x23),(x12*y21-x23*y32)*2]*ab/6; ...
                                [y12,0,x21];[0,x21,y12];[y12*(y32-y13),x21*(x23-x31),(x23*y32-x31*y13)*2]*ab/6]'/2/A;
    
    model(i).BL = BL;

    % Kb 
    Kb = A*BL'*Am*BL;
    
    % BNL

	Tx=1/(2*A)*[-y32 -y13 -y21];
	Ty=1/(2*A)*[-x23 -x31 -x12];
	
% 	T=[Tx;Ty];
    
    Bw=sparse(3,9);
	Bw(1,1)=1; Bw(2,4)=1; Bw(3,7)=1;
	
	Kxx=sparse(18,18); Kyy=sparse(18,18); Kxy=sparse(18,18);
	TBx = Tx*Bw;
    TBy = Ty*Bw;
    TByTBx = TBy'*TBx;
	Kxx([3 4 5 9 10 11 15 16 17],[3 4 5 9 10 11 15 16 17])=TBx'*TBx;
	Kyy([3 4 5 9 10 11 15 16 17],[3 4 5 9 10 11 15 16 17])=TBy'*TBy;
	Kxy([3 4 5 9 10 11 15 16 17],[3 4 5 9 10 11 15 16 17])= TByTBx + TByTBx';

    % contribution of (dv/dx)^2 and du/dy^2
    
%     Kxx([2 8 14],[2 8 14])=Tx'*Tx;
% 	Kyy([1 7 13],[1 7 13])=Ty'*Ty;

    % contribution of (du/dx)^2 and dv/dy^2
    
%     Kxx([1 7 13],[1 7 13])=Tx'*Tx;
% 	Kyy([2 8 14],[2 8 14])=Ty'*Ty;
    
    Bx = sparse(3,18); By = sparse(3,18);
    Bx(1,1)=1; Bx(2,7)=1; Bx(3,13)=1;
    By(1,2)=1; By(2,8)=1; By(3,14)=1;
    TxBx =  Tx*Bx;
    TyBx = Ty*Bx;
    TxBxTyBx = TxBx'*TyBx;
    TyBy = Ty*By;
    TxBy = Tx*By;
    TxByTyBy = TxBy'*TyBy;
    Kxy = Kxy + TxBxTyBx + TxBxTyBx' + TxByTyBy + TxByTyBy'; 
   
    model(i).Kxx = Kxx;
    model(i).Kyy = Kyy;
    model(i).Kxy = Kxy;
    
    
% if elements(i,5) == 1 
%     model(i).Kxx = 0*Kxx;
%     model(i).Kyy = 0*Kxy;
%     model(i).Kxy = 0*Kyy;
% end
% 
% if elements(i,5) == 3 
%     model(i).Kxx = 0*Kxx;
%     model(i).Kyy = 0*Kxy;
%     model(i).Kxy = 0*Kyy;
% end

    
    % higher order stiffness matrix
    
    % % hierarchical rotation to global displacement matrix
	Tqu=[[x32,y32,A4,x13,y13, 0,x21,y21, 0]; ...
         [x32,y32, 0,x13,y13,A4,x21,y21, 0]; ...
         [x32,y32, 0,x13,y13,0,x21,y21,A4]]/A4;
	
	% element sides length
	LL21=x21^2+y21^2; LL32=x32^2+y32^2; LL13=x13^2+y13^2;
	
	% cartesian to natural transformation matrix
	Te=[[y23*y13*LL21, y31*y21*LL32, y12*y32*LL13]; ...
        [x23*x13*LL21,x31*x21*LL32,x12*x32*LL13]; ...
        [(y23*x31+x32*y13)*LL21,(y31*x12+x13*y21)*LL32, (y12*x23+x21*y32)*LL13]]/(A*A4);
	
	% natural strain to hierarchical bending modes matrices
	Q1=[[b(1),b(2),b(3)]/LL21;[b(4),b(5),b(6)]/LL32;[b(7),b(8),b(9)]/LL13]*A2/3;
	Q2=[[b(9),b(7),b(8)]/LL21;[b(3),b(1),b(2)]/LL32;[b(6),b(4),b(5)]/LL13]*A2/3;
	Q3=[[b(5),b(6),b(4)]/LL21;[b(8),b(9),b(7)]/LL32;[b(2),b(3),b(1)]/LL13]*A2/3;
	Q4=(Q1+Q2)/2; Q5=(Q2+Q3)/2; Q6=(Q3+Q1)/2; 
	
	% natural stress-strain matrix
	Enat=Te'*Am*Te;
	
	% higher order stiffnes matrix in local coordinates (hierarchical bending modes)
	Kq=(3/4)*b0*A*(Q4'*Enat*Q4+Q5'*Enat*Q5+Q6'*Enat*Q6);
	
	% high order stiffness
	Kh=Tqu'*Kq*Tqu;
	
	% tangential stiffness matrix

	Kt([3 4 5 9 10 11 15 16 17],[3 4 5 9 10 11 15 16 17])=Allman_bending_stiffness(model(i).node,E,nu,h);
	% contribution for high order matrix
    Kt([1 2 6 7 8 12 13 14 18],[1 2 6 7 8 12 13 14 18])=Kh;
	
    model(i).Kt = Kt + Kb;
	   
end
