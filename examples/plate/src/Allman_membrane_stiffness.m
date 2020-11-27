function [K,Kb] = Allman_membrane_stiffness(nodes,E,nu,h)
% this function calculates the basic and higher-order stiffness matrix for a triangular mambrane element
% with three nodes and three dofs per node (2 displacement and 1 drilling rotation)
% 
% Ref: "A study of optimal membrane triangles with drilling freedoms"
%       C.A. Felippa
%       Report CU-CAS-03-02
%       Center for Aerospace Structures
%       University of Colorado
%
% x = x nodal coordinates [1x3]
% y = y nodal coordinates [1x3]
% h = thickness

% Am = material tensor (plain stress)
%
% A = element area
% K = element stiffness matrix (in global coordinates)

%nu = Am(1,2)/Am(1,1);

% set of parameters for the optimal membrane element
% ab = 3/2;
% b0 = (1-4*nu^2)/2;
% b = [1,2,1,0,1,-1,-1,-1,-2];

% for Allman's '88 membrane element (see Ref. page 14)
ab=1;
b0=4/9;
b = [1/12 5/12 1/2 0 1/3 -1/3 -1/12 -1/2 -5/12];

x=nodes(:,1); y=nodes(:,2);

Am=h*E/(1-nu^2)*[ 1 nu        0;
                  nu 1        0;
                  0  0 (1-nu)/2];
              
% geometric quantities
x12=x(1)-x(2); x23=x(2)-x(3); x31=x(3)-x(1); x21=-x12; x32=-x23; x13=-x31;
y12=y(1)-y(2); y23=y(2)-y(3); y31=y(3)-y(1); y21=-y12; y32=-y23; y13=-y31;

% element aera
A=(y21*x13-x21*y13)/2; A2=2*A; A4=4*A; 

% force lumping matrix
L= [[y23,0,x32];[0,x32,y23];[y23*(y13-y21),x32*(x31-x12),(x31*y13-x12*y21)*2]*ab/6; ...
    [y31,0,x13];[0,x13,y31];[y31*(y21-y32),x13*(x12-x23),(x12*y21-x23*y32)*2]*ab/6; ...
    [y12,0,x21];[0,x21,y12];[y12*(y32-y13),x21*(x23-x31),(x23*y32-x31*y13)*2]*ab/6]/2;


% basic stiffness matrix
Kb=(L*Am*L')/A;

% hierarchical rotation to global displacement matrix
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

% total stiffness matrix
K = Kh+Kb;

% numerical simmetry
K=1/2*(K'+K);