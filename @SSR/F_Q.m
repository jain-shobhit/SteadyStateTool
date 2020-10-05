function F_Q = F_Q(O,x)
% function to evaluate the zero function for quasiperiodic forcing
% S - function handle to the nonlinearity
switch O.sys_order
    case 'first'
        F = O.F_kappa - SSR_O1.evaluate_fun_over_array(O.R,z,false) * O.E;
        F_Q = O.Hmat*F(:) ;
    case 'second'
        
        F = O.f_kappa - SSR.evaluate_fun_over_array(O.S,x,false) * O.E;
        F_Q = O.Qmat*F(:) ;
end
