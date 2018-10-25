function update_f_kappa(O)
O.f_kappa = O.f_theta(O.theta_set)*O.E;
O.isupdated.f_kappa = true;
end