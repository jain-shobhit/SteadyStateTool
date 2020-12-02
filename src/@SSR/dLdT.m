function L = dLdT(O,t,T)
% t must be a row vector
% underdamped
om = repmat(O.omega(O.u), size(t));
al = repmat(O.alpha(O.u),size(t));
exp_alpha_T = repmat(exp(O.alpha(O.u) * T),size(t));
exp_alpha_2T = repmat(exp(2*O.alpha(O.u) * O.T),size(t));
exp_alpha_t_T_om = exp(O.alpha(O.u) * (t + T))./om;
denominator = repmat(1 + exp(2*O.alpha(O.u)*T) -...
    2 * exp(O.alpha(O.u)*T).* cos(O.omega(O.u)*T), size(t));
denominator2 = repmat((1 + exp(2*O.alpha(O.u)*T) -...
    2 * exp(O.alpha(O.u)*T).* cos(O.omega(O.u)*T)).^2, size(t));
Lu = exp_alpha_t_T_om.* ( ( om.* cos(O.omega(O.u)*(T+t)) + ...
    al.* sin(O.omega(O.u)*(T+t)) - 2*al.* exp_alpha_T.* sin(O.omega(O.u) * t) )./ ...
    denominator - ( (sin(O.omega(O.u)*(t + T)) - exp_alpha_T.* sin(O.omega(O.u) * t) ).* ...
    (2*(al.* (exp_alpha_2T - exp_alpha_T.* repmat(cos(O.omega(O.u)*T),size(t)) ) + ...
    om.* exp_alpha_T.* repmat(sin(O.omega(O.u)*T),size(t)) )) )./(denominator2));

% overdamped
if ~isempty(O.o)
    EaT = repmat(exp(O.gamma(O.o) * T), size(t));
    EatT = exp(O.gamma(O.o) * (t + T));
    EbT = repmat(exp(O.beta(O.o) * T), size(t));
    EbtT = exp(O.beta(O.o) * (t + T));
    EabT = repmat(exp((O.gamma(O.o) + O.beta(O.o)) * T), size(t));
    a = repmat(O.gamma(O.o) , size(t));
    b = repmat(O.beta(O.o) , size(t));
    denominator = (1 - EaT - EbT + EabT);
    Term1 = (b - (a + b).*EaT).* EbtT + ( (a + b).* EbT - a ).* EatT;
    Term2 =  ( EbtT.* (1 - EaT) - EatT.*(1 - EbT) ).*...
        ( a.*EaT + b.*EbT - (a + b).* EabT );
    Lo = ( Term1./denominator + Term2./(denominator.^2) )./(b-a);
else
    L = Lu;
    return
end
L = zeros(O.n,length(t));
L(O.u,:) = Lu;
L(O.o,:) = Lo;
end