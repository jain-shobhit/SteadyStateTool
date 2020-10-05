function [x, xd] = Picard(O,varargin)
% perform Picard iteration starting with optional initial guess
% optional maximum number of iterations, and optional tolerance

[x0,tol,maxiter] = parse_iteration_inputs(O,varargin);
%% Initialization
count = 1;
r = 1;
switch O.sys_order
    case 'first'
        z0 = x0; % fit notation with paper for first order implementation
        switch O.type
            case 'p'
                %% Iteration
                while r>tol
                    % Evaluate the map G_P
                    
                    R_array = SSR.evaluate_fun_over_array(O.R,z0,false);
                    eta = O.convolution_x( O.Vinv * (O.Ft - R_array) );
                    z = O.V * eta;
                    % convergence check
                    r = norm(z-z0,Inf)/norm(z0,Inf);
                    %disp(['Iteration ' num2str(count) ': ' 'relative change = ' num2str(r)])
                    if count>maxiter || isnan(r)
                        error('Picard iteration did not converge')
                    end
                    z0 = z;
                    count = count + 1;
                end
            case 'qp'
                while r>tol
                    R_array = SSR.evaluate_fun_over_array(O.R,z0,false);
                    F = O.F_kappa - R_array*O.E;
                    z_kappa = reshape(O.Hmat * F(:),O.N,[]);
                    z = real(z_kappa*O.Einv);
                    % convergence check
                    r = norm(z-z0,Inf)/norm(z0,Inf);
                    %disp(['Iteration ' num2str(count) ': ' 'relative change = ' num2str(r)])
                    if count>maxiter || isnan(r)
                        error('Picard iteration did not converge')
                    end
                    z0 = z;
                    count = count + 1;
                end
        end
        x = z;
        xd = []; %velocities are already stored in z
        
    case 'second'
        switch O.type
            case 'p'
                %% Iteration
                while r>tol
                    % Evaluate the map G_P
                    S_array = SSR.evaluate_fun_over_array(O.S,x0,false);
                    eta = O.convolution_x( O.U.' * (O.Ft - S_array) );
                    x = O.U * eta;
                    % convergence check
                    r = norm(x-x0,Inf)/norm(x0,Inf);
                    disp(['Iteration ' num2str(count) ': ' 'relative change = ' num2str(r)])
                    if count>maxiter || isnan(r)
                        error('Picard iteration did not converge')
                    end
                    x0 = x;
                    count = count + 1;
                end
                xd = O.U * O.convolution_xd( O.U.' * (O.Ft - S_array) );
            case 'qp'
                while r>tol
                    S_array = SSR.evaluate_fun_over_array(O.S,x0,false);
                    F = O.f_kappa - S_array*O.E;
                    x_kappa = reshape(O.Qmat * F(:),O.n,[]);
                    x = real(x_kappa*O.Einv);
                    % convergence check
                    r = norm(x-x0,Inf)/norm(x0,Inf);
                    disp(['Iteration ' num2str(count) ': ' 'relative change = ' num2str(r)])
                    if count>maxiter || isnan(r)
                        error('Picard iteration did not converge')
                    end
                    x0 = x;
                    count = count + 1;
                end
                xd = real(x_kappa*O.Evinv);
        end
end