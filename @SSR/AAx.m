function AA_zeta = AAx(O,x0)
% compute the integral Ax via the precomputed integrals in Lint
Lint = get_Lint(O);
AA_zeta = squeeze(sum(Lint.*x0,2));
end