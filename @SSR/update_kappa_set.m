function update_kappa_set(O)
n_freq = length(O.Omega);
% check inputs
if length(O.nh) ~= n_freq
    if length(O.nh) == 1
        nmax = repmat(O.nh,n_freq,1);
    else
        error(['dimenionality of number of harmonics entered' ...
            ' does not match the dimension of base frequency. '...
            'Please check input dimensions']);
    end
else
    nmax = O.nh;
end
% compute kappas
range = cell(n_freq,1);
for j = 1:n_freq
    range{j} = -nmax(j):nmax(j);
end
O.kappa_set = combvec(range{:});
O.isupdated.kappa = true;
end