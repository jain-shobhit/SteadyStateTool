function J = J(O,t,T)
% t must be a row vector
J = zeros(O.n,length(t));
% underdamped
om = repmat(O.omega(O.u), size(t));
al = repmat(O.alpha(O.u), size(t));

exp_alpha_T = repmat(exp(O.alpha(O.u) * T),size(t));
exp_alpha_t_om = exp(O.alpha(O.u) * t)./om;
exp_alpha_t_T_om = exp(O.alpha(O.u) * (t + T))./om;
hu = repmat(SSR.h(t),size(O.omega(O.u)));
denominator = repmat(1 + exp(2*O.alpha(O.u)*T) -...
    2 * exp(O.alpha(O.u)*T).* cos(O.omega(O.u)*T), size(t));
J(O.u,:) = (exp_alpha_t_T_om.* ( om.*(cos(O.omega(O.u)*(t+T)) - ...
    exp_alpha_T.* cos(O.omega(O.u)*t)) +...
    al.*(sin(O.omega(O.u)*(O.T+t)) - ...
    exp_alpha_T.* sin(O.omega(O.u)*t))))./denominator ...
    + hu.* (exp_alpha_t_om.* (om.*cos(O.omega(O.u)*t) + ...
    al.*sin(O.omega(O.u)*t)));

% critically damped
if ~isempty(O.c)
    t_c = repmat(t,size(O.c));
    al = repmat(O.alpha(O.c), size(t));
    exp_alpha_t = exp(O.alpha(O.c) * t);
    exp_alpha_t_T = exp(O.alpha(O.c) * (t + T));
    one_min_exp_alpha_T = 1 - repmat(exp(O.alpha(O.c) * T),size(t));
    hc = repmat(SSR.h(t),size(O.omega(O.c)));
    J(O.c,:) = exp_alpha_t_T.*(one_min_exp_alpha_T.*(1 + al.*t_c) +...
        al.*T)./(one_min_exp_alpha_T.^2) + hc.*exp_alpha_t.*(1 + ...
        al.*t_c) ;
end

% overdamped
if ~isempty(O.o)
    ho = repmat(SSR.h(t),size(O.omega(O.o)));
    beta_min_gamma = repmat(O.beta(O.o) - O.gamma(O.o),size(t));
    bet = repmat(O.beta(O.o),size(t));
    gam = repmat(O.gamma(O.o),size(t));
    exp_beta_t_T = exp(O.beta(O.o) * (t + T));
    exp_gamma_t_T = exp(O.gamma(O.o) * (t + T));
    exp_beta_T = repmat(exp(O.beta(O.o) * T),size(t));
    exp_gamma_T = repmat(exp(O.gamma(O.o) * T),size(t));
    exp_beta_t = exp(O.beta(O.o) * t);
    exp_gamma_t = exp(O.gamma(O.o) * t);
    
    J(O.o,:) = ( (bet.*exp_beta_t_T./(1 - exp_beta_T) - ...
        gam.*exp_gamma_t_T./(1 - exp_gamma_T)) + ...
        ho.*(bet.*exp_beta_t - gam.*exp_gamma_t) )./ beta_min_gamma;
end
end