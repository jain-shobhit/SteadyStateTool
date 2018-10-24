function [x, xd] = NewtonRaphson(O,varargin)
% perform Newton--Raphson iteration starting with an optional
% initial guess, optional maximum number of iterations and
% optional tolerance
[x0,tol,maxiter] = parse_iteration_inputs(O,varargin);

count = 1;
r=1;

switch O.domain
    case 'time'
        while r>tol
            F = F_P(O,x0);      % evaluation of the zero function
            DF = DF_P(O,x0);    % derivative of zero function
            mu = -DF\F(:);      % correction to the solution guess
            x = x0 + reshape(mu,size(x0));  % new solution guess
            r = norm(x-x0,Inf)/norm(x0,Inf);
            disp(['Iteration ' num2str(count) ': ' '||dx||/||x|| = ' num2str(r), ' residual = ' num2str(norm(F(:)))])
            if count>maxiter || isnan(r)
                warning('Newton iterations did not converge')
                x = [];
                xd = [];
                return
            end
            x0 = x;
            count = count + 1;
        end
        S_array = SSR.evaluate_fun_over_array(O.S,x,false);
        xd = O.U * O.convolution_xd(O.U.'*(O.Ft - S_array));
        
    case 'freq'
        x_kappa0 = x0*O.E;                          % Fourier domain
        while r>tol
            F = x_kappa0(:) - F_Q(O,x0);                % evaluation of the zero function
            DF = DF_Q(O,x0);                            % derivative of zero function
            mu = -DF\F(:);                              % correction to the solution guess
            x_kappa = x_kappa0 + reshape(mu,O.n,[]);    % new solution guess
            x = real(x_kappa * O.Einv);                 % solution in time domain
            r = norm(x - x0,Inf)/norm(x0,Inf);
            disp(['Iteration ' num2str(count) ': ' '||dx||/||x|| = ' num2str(r), ' residual = ' num2str(norm(F(:)))])
            if count>maxiter || isnan(r)
                warning('Newton iterations did not converge')
                x = [];
                xd = [];
                return
            end
            x_kappa0 = x_kappa;
            x0 = x;
            count = count + 1;
        end
        xd = real(x_kappa*O.Evinv);
end
end