function [Mblu,mass]=Allman_bending_mass(nodes,h,rho)

% Allman triangle bending consistent mass matrix
% 
% Ref: D.J. Allman
%"Implementation of a flat facet shell finite element for application in structural dynamics"
% Computers & Structures, vol 59, pp 657-663

co1=nodes(1,:);
co2=nodes(2,:);
co3=nodes(3,:);

x1=co1(1); y1=co1(2);
x2=co2(1); y2=co2(2);
x3=co3(1); y3=co3(2);

x12=x1-x2;  y21=y2-y1;
x23=x2-x3;  y32=y3-y2;
x31=x3-x1;  y13=y1-y3;

% N=[e1 e2 e3];
% Nalpha=[e1*e2 e2*e3 e3*e1 e1*e2*(e2-e1) e2*e3*(e3-e2) e3*e1*(e1-e3)];

Bw=zeros(3,9);
Bw(1,1)=1; Bw(2,4)=1; Bw(3,7)=1;

Baw=zeros(6,9);
Baw(1,[2 3 5 6])=1/2*[y21 x12 -y21 -x12];
Baw(2,[5 6 8 9])=1/2*[y32 x23 -y32 -x23];
Baw(3,[2 3 8 9])=1/2*[-y13 -x31 y13 x31];
Baw(4,1:6)    =1/2*[-2 -y21 -x12 2 -y21 -x12];
Baw(5,4:9)    =1/2*[-2 -y32 -x23 2 -y32 -x23];
Baw(6,[1:3 7:9])=1/2*[2 -y13 -x31 -2 -y13 -x31];

co=[1  1  1;
    x1 x2 x3;
    y1 y2 y3];
    
A=1/2*det(co);

intNN=[1/6 1/12 1/12;
       0   1/6   1/12;
       0     0   1/6];

% intNN=rho*h*A*1/2*(intNN'+intNN);
intNN=rho*h*A*(intNN+tril(intNN',-1));

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

 %intNaNa=1/2*(intNaNa+intNaNa');
 intNaNa=(intNaNa+tril(intNaNa',-1));
 
 Mb=Bw'*intNN*Bw+Bw'*intNNa*Baw+Baw'*intNaN*Bw+Baw'*intNaNa*Baw;
 
 mass = A*h*rho;

s = Mb(1,1) + Mb(4,4) + Mb(7,7);

lumped = diag(Mb).*(mass./s);

Mblu = diag(lumped);