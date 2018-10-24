function update_Qmat(O)
kappas = O.kappa_set;
n_kappa = size(kappas,2);
Q_T_kappa = zeros(O.n,n_kappa);

% underdamped modes
Q_T_kappa(O.u,:) = 1./( (1i* (O.Omega*kappas) ...
    - repmat(O.alpha(O.u),1,n_kappa)).^2 + ...
    repmat(O.omega(O.u).^2,1,n_kappa) );
% critically damped modes
if ~isempty(O.c)
    Q_T_kappa(O.c,:) = 1./(1i * (O.Omega*kappas) ...
        - repmat(O.alpha(O.c),1,n_kappa) ).^2;
end
% overdamped modes
if ~isempty(O.o)
    Q_T_kappa(O.o,:) = 1./( (repmat(O.beta(O.o),1,n_kappa) - ...
        1i * (O.Omega*kappas)).*(repmat(O.gamma(O.o),1,n_kappa) - ...
        1i * (O.Omega*kappas)) );
end

%             N = O.n*n_kappa; % size of the matrix
%             Ublk = kron(speye(n_kappa,n_kappa),O.U);
%             O.Qmat = Ublk*spdiags(Q_T_kappa(:),0,N,N)*Ublk.';

Q_kappa = cell(n_kappa,1);
for j = 1:n_kappa
    Q = Q_T_kappa(:,j);
    Q_kappa{j} = O.U * diag(Q) * O.U.';
end
O.Qmat = spblkdiag(Q_kappa{:});
O.isupdated.Q = true;
end