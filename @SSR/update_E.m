function update_E(O)
kappa = O.kappa_set;
theta = O.theta_set;
n_theta = size(theta,2);
N = size(theta,2); % total number of elements on the torus
thetakappa = 2i*pi*theta.'*kappa;
EXP = exp(-thetakappa)/N;
O.E = EXP;
O.Ekron = kron(EXP.',speye(O.n, O.n));
EXPinv = exp(thetakappa.');
O.Einv = EXPinv;
O.Evinv = EXPinv.*repmat(1i*(O.Omega*kappa).',1,n_theta);
O.Einvkron = kron(EXPinv.',speye(O.n, O.n));
O.isupdated.E = true;
end