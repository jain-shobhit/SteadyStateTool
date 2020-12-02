function [y, yd, z0] = NewtonRaphson_reformulated(O,varargin)
% perform Newton--Raphson iteration starting with an optional
% initial guess, optional maximum number of iterations and
% optional tolerance
[z0,eta0,tol,maxiter] = parse_iteration_inputs_ref(O,varargin);
count = 1;
r = 1;
switch O.type
    case 'p'
        while r>tol
            [F,Tz0] = F_P_reformulated(O,eta0,z0);      % evaluation of the zero function
            UTz0 = O.U * Tz0;
            DF = DF_P_reformulated(O,UTz0);             % derivative of zero function
            
            % check if the problem is already parallelized by DF_P
            if isdistributed(DF)
                F2 = distributed(F(:));
                spmd
                    mu = -DF\F2;      % correction to the solution guess
                end
                mu = gather(mu);
            else
                mu = -DF\F(:);        % correction to the solution guess
            end

            z = z0 + reshape(mu,size(z0));  % new solution guess
            Tz = O.AAx(z) + eta0;
            UTz = O.U * Tz;
            r = norm(UTz - UTz0,Inf)/norm(UTz0,Inf);
            disp(['Iteration ' num2str(count) ': ' '||dx||/||x|| = ' num2str(r), ' residual = ' num2str(norm(F(:)))])
            if count>maxiter || isnan(r)
                warning('Newton iterations did not converge')
                y = [];
                yd = [];
                return
            end
            z0 = z;
            UTz0 = UTz;
            count = count + 1;
        end 
        y = UTz0;
        S_array = SSR.evaluate_fun_over_array(O.S,UTz0,false);
        yd = O.U * O.convolution_xd(O.U.'*(O.Ft - S_array));
    case 'qp'
        error('Reformulation is not implemented for quasiperiodic systems!')
end
end