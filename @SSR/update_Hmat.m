function update_Hmat(O)

kappas = O.kappa_set;
n_kappa = size(kappas,2);

H_T_kappa = 1./( 1i* O.Omega*kappas ...
    - repmat(O.lambda,1,n_kappa)); % size: O.N, n_kappa


H_kappa = cell(n_kappa,1);
for j = 1:n_kappa
    H = H_T_kappa(:,j);
    H_kappa{j} = O.V * diag(H) * O.Vinv;
end

if verLessThan('matlab','9.7.0')
    O.Hmat = spblkdiag(H_kappa{:});
else
    O.Hmat = blkdiag(H_kappa{:});
end 

O.isupdated.H = true;
end