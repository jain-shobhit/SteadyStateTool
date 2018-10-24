function F_P = F_P(O,x)
% function to evaluate the zero function for periodic forcing
% eta  - modal unknowns
% S - function handle to the nonlinearity
% f - external forcing
S_array = SSR.evaluate_fun_over_array(O.S,x,false);
F_P = x + O.U * O.convolution_x(O.U.'*(S_array-O.Ft));
end