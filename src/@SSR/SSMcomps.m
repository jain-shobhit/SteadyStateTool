function SSMcomps(O)
% computes additional material required for SSM computation
[V1, dd] = eig(full(O.K),full(O.M),'vector');
[~, ind] = sort(dd);
V1 = V1(:,ind);
V2 = V1(:,setdiff(1:length(O.M),O.mode_choice));
mu2 = diag(V2.' * O.M * V2);
U2 = V2 * diag( 1./ sqrt(mu2) );
omega2 = sqrt(diag((U2.' * O.K * U2)));
zeta2 = diag((U2.' * O.C * U2))./ (2*omega2);
O.U2store = U2;
O.omega2store = omega2;
O.zeta2store = zeta2;
O.SSMdone = 1;
end