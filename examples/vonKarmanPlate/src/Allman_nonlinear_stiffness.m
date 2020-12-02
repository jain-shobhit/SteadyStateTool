function [K2,K3] = Allman_nonlinear_stiffness(BL,Kxx,Kyy,Kxy,~,Am,A)
% Am - 3X3
% A - S
% Kt - 18X18
% qe - 18X1
% Kxy - 18X18
% Kxx - 18X18
% Kyy - 18X18
% BNL - 3X18
% N - 3X1
% BL - 3X18
% Note - symmetry of Kxx, Kyy, Kxy used while calculating non linear
% tensors, Symmetry of Am not taken into account.

nndof = 18 ;
% F = Kt*qe;
% % N
% BNL=[qe'*Kxx;
%     qe'*Kyy;
%     qe'*Kxy];
% 
% N=Am*(BL+1/2*BNL)*qe;
% 
% % Kt=Kt+A*(BL+BNL)'*Am*(BL+BNL)+A*(N(1)*Kxx+N(2)*Kyy+N(3)*Kxy);
% 
% %Kt=A*(BL+BNL)'*Am*(BL+BNL)+A*(N(1)*Kxx+N(2)*Kyy+N(3)*Kxy);
% 
% Kt=Kt+A*(BL'*Am*BNL + BNL'*Am*BL + BNL'*Am*BNL)+A*(N(1)*Kxx + N(2)*Kyy + N(3)*Kxy);
% % Internal forces
% F= F + A*(BL+BNL)'*N - A*BL'*Am*BL*qe;

%% Non linear element stiffness tensors calculation

K2 = zeros(nndof,nndof,nndof);
K3 = zeros(nndof,nndof,nndof,nndof);

C = A*BL'*Am;
D = A*Am*BL; % D = C'
for i = 1:nndof
    K2_1 = 1/2*(C(i,1)*(Kxx) + C(i,2)*(Kyy) + C(i,3)*(Kxy)) ;             % corresponds to A*BL'*Am*BNL
    K2_2 = Kxx(:,i)*D(1,:) + Kyy(:,i)*D(2,:) + Kxy(:,i)*D(3,:);     % corresponds to BNL'*Am*BL*A
    K2(i,:,:) = K2_1 + K2_2 ;
    for j = 1:nndof
        K3_1 = 1/2*A* (( Kxx(i,j)*Am(1,1)+ Kyy(i,j)*Am(2,1)+Kxy(i,j)*Am(3,1) )*Kxx + ...   % corresponds to A*BNL'*Am*BNL
                        ( Kxx(i,j)*Am(1,2)+ Kyy(i,j)*Am(2,2)+Kxy(i,j)*Am(3,2) )*Kyy + ...
                        ( Kxx(i,j)*Am(1,3)+ Kyy(i,j)*Am(2,3)+Kxy(i,j)*Am(3,3) )*Kxy) ;
        
        K3(i,j,:,:) = K3_1;
    end
end



