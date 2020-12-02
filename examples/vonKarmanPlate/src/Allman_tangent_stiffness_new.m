function [F,Kt] = Allman_tangent_stiffness_new(BL,Kxx,Kyy,Kxy,Kt,Am,A,qe)

F = Kt*qe;
% N
BNL=[qe'*Kxx;
     qe'*Kyy;
     qe'*Kxy];

N=Am*(BL+1/2*BNL)*qe;
Kt=Kt+A*(BL'*Am*BNL+BNL'*Am*BL+BNL'*Am*BNL)+A*(N(1)*Kxx+N(2)*Kyy+N(3)*Kxy);
% Internal forces
F= F + A*(BL+BNL)'*N - A*BL'*Am*BL*qe;

