function [x,xd] = LinearResponse(O)

switch O.sys_order
    case 'first'
        switch O.type
            case 'p'
                % Compute modal forcing
                phi = O.W' * O.Ft;
                % Compute modal response by convolution with Green's function
                eta = convolution_x(O,phi);
                z = O.V*eta;          
            case 'qp'
                
                z_kappa = reshape(O.Hmat * O.F_kappa(:),O.N,[]);
                z = real(z_kappa*O.Einv);      
        end
        x  = z; %Outputs the vector z that contains positions and velocities
        xd = [];
    case 'second'
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