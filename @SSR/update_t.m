function update_t(O)
if O.order
    O.nt = O.n_steps*O.order + 1; % first set n_steps
    O.dt = O.T/(O.nt - 1);        % time node spacing
    w = SSR.NewtonCotes(O.order); % weights for integration over an interval
    w0 = [w(2:end-1) w(1)+w(end)]; % repeating block in the weights vector
    O.weights = O.dt*[w(1) repmat(w0,1,O.n_steps-1) w(2:end)]; % weights for the whole grid
else
    O.nt = O.n_steps + 1;
    O.dt = O.T/(O.nt - 1);
    O.weights = O.dt;
end

O.t = 0:O.dt:O.T;
not_updated(O);
end