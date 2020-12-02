function [w,VN] = sortVMs(V,D)

% Description:
%--------------------------------------------------------------------------
% Order the eigenvectors and corresponding frequencies obtained from
% solving the eigenvalue problem.
%
% Input:
% VMs       : Eigenvectors
% D         : Eigenvalues (omega^2)

% Output:
% f         : Eigenfrequency [Herz]
% w         : Eigenfrequency [rad/s]
% VN        : Eigevectors

%--------------------------------------------------------------------------

% This function checks if the eigenvalues and eigenvectors are ordered
% starting with the lowest eigenvalue solution, if not it will sort it to
% this form.

w = (diag(D));
[w,pos] = sort(w);
VN = V(:,pos);

end
