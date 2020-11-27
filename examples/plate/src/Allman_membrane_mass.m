function [Mmlu,mass]=Allman_membrane_mass(nodes,h,rho)

% Allman triangle membrane consistent mass matrix
% 
% Ref: D.J. Allman
%"Implementation of a flat facet shell finite element for applcication in structural dynamics"
% Computers & Structures, val 59, pp 657-663

%syms e1 e2 e3 real

co1=nodes(1,:);
co2=nodes(2,:);
co3=nodes(3,:);

x1=co1(1); y1=co1(2);
x2=co2(1); y2=co2(2);
x3=co3(1); y3=co3(2);

x12=x1-x2;  y21=y2-y1;
x23=x2-x3;  y32=y3-y2;
x31=x3-x1;  y13=y1-y3;

co=[1  1  1;
    x1 x2 x3;
    y1 y2 y3];
    
A=1/2*det(co);

X12=-(x1-x2)/(4*A);  Y12=-(y1-y2)/(4*A);
X23=-(x2-x3)/(4*A);  Y23=-(y2-y3)/(4*A);
X31=-(x3-x1)/(4*A);  Y31=-(y3-y1)/(4*A);

% N=[e1 e2 e3];
% Nalpha=[e1*e2 e2*e3 e3*e1 e1*e2*(e2-e1) e2*e3*(e3-e2) e3*e1*(e1-e3)];

Bu=zeros(3,9);
Bu(1,1)=1; Bu(2,4)=1; Bu(3,7)=1;

Bv=zeros(3,9);
Bv(1,2)=1; Bv(2,5)=1; Bv(3,8)=1;

Bau=zeros(6,9);

Bau(1,[3 6])=1/2*[-y21 y21];
Bau(2,[6 9])=1/2*[-y32 y32];
Bau(3,[3 9])=1/2*[y13 -y13];
Bau(4,1:8)=[X23*y21 Y23*y21 1/2*y21 X31*y21 Y31*y21 1/2*y21 X12*y21 Y12*y21];
Bau(5,[1:2 4:9])=[X23*y32 Y23*y32 X31*y32 Y31*y32 1/2*y32 X12*y32 Y12*y32 1/2*y32];
Bau(6,[1:5 7:9])=[X23*y13 Y23*y13 1/2*y13 X31*y13 Y31*y13 X12*y13 Y12*y13 1/2*y13];


Bav=zeros(6,9);

Bav(1,[3 6])=1/2*[-x12 x12];
Bav(2,[6 9])=1/2*[-x23 x23];
Bav(3,[3 9])=1/2*[x31 -x31];
Bav(4,1:8)=[X23*x12 Y23*x12 1/2*x12 X31*x12 Y31*x12 1/2*x12 X12*x12 Y12*x12];
Bav(5,[1:2 4:9])=[X23*x23 Y23*x23 X31*x23 Y31*x23 1/2*x23 X12*x23 Y12*x23 1/2*x23];
Bav(6,[1:5 7:9])=[X23*x31 Y23*x31 1/2*x31 X31*x31 Y31*x31 X12*x31 Y12*x31 1/2*x31];


intNN=[1/6 1/12 1/12;
       0   1/6   1/12;
       0     0   1/6];

intNN=rho*h*A*1/2*(intNN'+intNN);

intNNa=rho*h*A*[1/30 1/60 1/30 -1/180      0  -1/180;
                1/30 1/30 1/60  1/180 -1/180       0;
                1/60 1/30 1/30      0  1/180  -1/180];
            
intNaN=intNNa';

intNaNa=h*rho*A*[1/90 1/180 1/180       0 -1/1260  1/1260;
                    0  1/90 1/180  1/1260       0 -1/1260;
                    0     0  1/90 -1/1260  1/1260       0;
                    0     0     0   1/840 -1/2520 -1/2520;
                    0     0     0       0   1/840 -1/2520;
                    0     0     0       0       0   1/840];

 intNaNa=1/2*(intNaNa+intNaNa');

% Mu=rho*h*(Bu'*intNN*Bu+Bu'*intNNa*Bau+Bau'*intNaN*Bu+Bau'*intNaNa*Bau);
Mu=(Bu'*intNN*Bu+Bu'*intNNa*Bau+Bau'*intNaN*Bu+Bau'*intNaNa*Bau);


% Mv=rho*h*(Bv'*intNN*Bv+Bv'*intNNa*Bav+Bav'*intNaN*Bv+Bav'*intNaNa*Bav);
Mv=(Bv'*intNN*Bv+Bv'*intNNa*Bav+Bav'*intNaN*Bv+Bav'*intNaNa*Bav);
Mm=Mu+Mv;

mass = A*h*rho;

s = Mm(2,2) + Mm(5,5) + Mm(8,8);

lumped = diag(Mm).*(mass./s);

Mmlu = diag(lumped);


