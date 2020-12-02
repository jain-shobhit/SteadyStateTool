function [F_P,Tz] = F_P_reformulated(O,eta0,z)
% function to evaluate the zero function for periodic forcing
% eta0  - linear response; same as f_m
% S - function handle to the nonlinearity
A_zeta = O.AAx(z);               % Az
Tz =  (eta0 + A_zeta);          % Tz 
S_array = SSR.evaluate_fun_over_array(O.S,O.U * Tz,false);
F_P = z + O.U.' * S_array;      % the zero function z + GTz
end