function G = G(O,t,T)
% Compute greens function as in (9)
% t must be a row vector
% O.lambda is a column vector with all eigenvalues

exp_lambda_t = exp(O.lambda *t);
exp_lambda_T = exp(O.lambda*T) ;
fraction     = exp_lambda_T ./ ( 1- exp_lambda_T);
G = exp_lambda_t .*( fraction + SSR.h(t)) ;

end