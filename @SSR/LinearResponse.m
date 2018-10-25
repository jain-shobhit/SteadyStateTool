function [x,xd] = LinearResponse(O)
switch O.type
    case 'p'
        % Compute modal forcing
        phi = O.U.' * O.Ft;
        % Compute modal response by convolution with Green's function
        eta = convolution_x(O,phi);
        etad = convolution_xd(O,phi);
        x = O.U*eta;
        xd = O.U*etad;
    case 'qp'
        x_kappa = reshape(O.Qmat * O.f_kappa(:),O.n,[]);
        x = real(x_kappa*O.Einv);
        xd = real(x_kappa*O.Evinv);
end
end