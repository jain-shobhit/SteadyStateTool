function F_P = F_P(O,x)
% function to evaluate the zero function for periodic forcing
% eta  - modal unknowns
% S - function handle to the nonlinearity
% f - external forcing

switch O.sys_order
    case 'first'
        z = x;
        R_array = SSR.evaluate_fun_over_array(O.R,z,false);
        F_P = z + O.V * O.convolution_x(O.W' *(R_array-O.Ft));
    case 'second'
        B = size(O.U);
        if B(1) == B(2)
            S_array = SSR.evaluate_fun_over_array(O.S,x,false);
            F_P = x + O.U * O.convolution_x(O.U.'*(S_array-O.Ft));
        else
            S_array = SSR.evaluate_fun_over_array(O.S,O.U * x,false);
            F_P = x + O.convolution_x(O.U.'*(S_array-O.Ft));
        end
end