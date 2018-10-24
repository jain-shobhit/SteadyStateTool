function [x, xd] = Picard(O,varargin)
% perform Picard iteration starting with optional initial guess
% optional maximum number of iterations, and optional tolerance

[x0,tol,maxiter] = parse_iteration_inputs(O,varargin);
%% Initialization
count = 1;
r = 1;
switch O.domain
    case 'time'
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
    case 'freq'
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