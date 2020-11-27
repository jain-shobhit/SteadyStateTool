function [Kt] = Allman_tangent_stiffness_derv(BL,Kxx,Kyy,Kxy,~,Am,A,qe)
% Evaluating sensitivity at Equilibrium position ( Only linear terms )
BNL=[qe'*Kxx;
     qe'*Kyy;
     qe'*Kxy];

% N=Am*(BL+1/2*BNL)*qe;
N=Am*BL*qe;

% Kt=Kt+A*(BL'*Am*BNL+BNL'*Am*BL+BNL'*Am*BNL)+A*(N(1)*Kxx+N(2)*Kyy+N(3)*Kxy);
Kt = A*(BL'*Am*BNL + BNL'*Am*BL )+A*(N(1)*Kxx + N(2)*Kyy + N(3)*Kxy);

