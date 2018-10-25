function F_Q = F_Q(O,x)
% function to evaluate the zero function for quasiperiodic forcing
% S - function handle to the nonlinearity
F = O.f_kappa - SSR.evaluate_fun_over_array(O.S,x,false) * O.E;
F_Q = O.Qmat*F(:) ;
end
