function update_f_kappa(O)

switch O.sys_order
    case 'first'
        O.F_kappa = O.F_theta(O.theta_set)*O.E;
        O.isupdated.F_kappa = true;
    case 'second'
        
        O.f_kappa = O.f_theta(O.theta_set)*O.E;
        O.isupdated.f_kappa = true;
end