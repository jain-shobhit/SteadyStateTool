function [E] = GRE(u,u1,M,K)

% R = diag(sqrt((u - u1)'*(u - u1)))/sqrt(max(diag(u'*u)));      
E.R = 'Not Calculated'; 
ux = u(1:6:end,:);
ux1 = u1(1:6:end,:);
E.GREx = norm(ux-ux1,'fro')/norm(ux,'fro');
% GREx = sqrt(sum(diag((u(1:6:end,:) - u1(1:6:end,:))'*(u(1:6:end,:) - ...
%     u1(1:6:end,:))))) / (sqrt(sum(diag(u(1:6:end,:)'*u(1:6:end,:)))));
uy = u(2:6:end,:);
uy1 = u1(2:6:end,:);
E.GREy = norm(uy-uy1,'fro')/norm(uy,'fro');

% GREy = sqrt(sum(diag((u(2:6:end,:) - u1(2:6:end,:))'*(u(2:6:end,:) - ...
%     u1(2:6:end,:))))) / (sqrt(sum(diag(u(2:6:end,:)'*u(2:6:end,:)))));
uz = u(3:6:end,:);
uz1 = u1(3:6:end,:);
E.GREz = norm(uz-uz1,'fro')/norm(uz,'fro');
 
% GREz = sqrt(sum(diag((u(3:6:end,:) - u1(3:6:end,:))'*(u(3:6:end,:) - ...
%     u1(3:6:end,:))))) / (sqrt(sum(diag(u(3:6:end,:)'*u(3:6:end,:)))));
Merr = M*(u-u1);
uu1 = u-u1;
Mu = M*u;
E.GRE_M = sqrt(uu1(:).' * Merr(:)) / sqrt(u(:).' * Mu(:)); 

Kerr = K*(u-u1);
Ku = K*u;
E.GRE_K = sqrt(uu1(:).' * Kerr(:)) / sqrt(u(:).' * Ku(:)); 