function [x, xd] = NewtonRaphson(O,varargin)
% perform Newton--Raphson iteration starting with an optional
% initial guess, optional maximum number of iterations and
% optional tolerance
[x0,tol,maxiter] = parse_iteration_inputs(O,varargin);

count = 1;
r=1;
switch O.sys_order
    case 'first'
        z0 = x0; % explicit assignment for clearness
        
        switch O.type
            case 'p'
                while r>tol
                    F = F_P(O,z0);      % evaluation of the zero function
                    DF = DF_P(O,z0);    % derivative of zero function
                    mu = -lsqminnorm(DF,F(:));      % correction to the solution guess
                    z = z0 + reshape(mu,size(z0));  % new solution guess
                    r = norm(z-z0,Inf)/norm(z0,Inf);
                    %disp(['Iteration ' num2str(count) ': ' '||dx||/||x|| = ' num2str(r), ' residual = ' num2str(norm(F(:)))])
                    if count>maxiter || isnan(r)
                        warning('Newton iterations did not converge')
                        x = [];
                        xd = [];
                        return
                    end
                    z0 = z;
                    count = count + 1;
                end
                
            case 'qp'
                z_kappa0 = z0*O.E;                          % Fourier domain
                while r>tol
                    F = z_kappa0(:) - F_Q(O,z0);                % evaluation of the zero function
                    DF = DF_Q(O,z0);                            % derivative of zero function
                    mu = -lsqminnorm(DF,F(:));                  % correction to the solution guess
                    
                    z_kappa = z_kappa0 + reshape(mu,O.N,[]);    % new solution guess
                    z = real(z_kappa * O.Einv);                 % solution in time domain
                    r = norm(z - z0,Inf)/norm(z0,Inf);
                    %disp(['Iteration ' num2str(count) ': ' '||dx||/||x|| = ' num2str(r), ' residual = ' num2str(norm(F(:)))])
                    if count>maxiter || isnan(r)
                        warning('Newton iterations did not converge')
                        x = [];
                        xd = [];
                        return
                    end
                    z_kappa0 = z_kappa;
                    z0 = z;
                    count = count + 1;
                end
        end
        
        x  = z; % Explicit assignment for clearness
        xd = []; % velocities are contained in z
        
        
    case 'second'
        
        switch O.type
            case 'p'
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
                
            case 'qp'
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