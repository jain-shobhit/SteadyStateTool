function update_theta_set(O)
n_freq = length(O.Omega);
% check inputs
if length(O.m) ~= n_freq
    if length(O.m) == 1
        mmax = repmat(O.m,n_freq,1);
    else
        error(['dimenionality of max. grid points number entered' ...
            ' does not match the dimension of base frequency. '...
            'Please check input dimensions']);
    end
else
    mmax = O.m;
end
% compute kappas
range = cell(n_freq,1);
for j = 1:n_freq
    dtheta = 1/mmax(j);
    range{j} = 0:dtheta:(1-dtheta);
end
O.theta_set = combvec(range{:});
O.isupdated.theta = true;
end