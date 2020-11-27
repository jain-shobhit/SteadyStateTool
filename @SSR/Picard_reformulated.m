function [y, yd, z0] = Picard_reformulated(O,varargin)
% perform Picard iteration starting with optional initial guess
% optional maximum number of iterations, and optional tolerance

[z0,eta0,tol,maxiter] = parse_iteration_inputs_ref(O,varargin);
%% Initialization
count = 1;
r = 1;
switch O.type
    case 'p'
        %% Iteration 
        Ueta0 = O.U * eta0;
        while r>tol
            % Evaluate the map GT
            A_zeta = O.AAx(z0);
            Ueta = O.U * (A_zeta + eta0);
            S_eta = SSR.evaluate_fun_over_array(O.S,Ueta,false);
            z = - O.U.' * S_eta;
            % convergence check
            r = norm(Ueta - Ueta0,Inf)/norm(Ueta0,Inf);
            disp(['Iteration ' num2str(count) ': ' 'relative change = ' num2str(r)])
            if count>maxiter || isnan(r)
                error('Picard iteration did not converge')
            end
            z0 = z;
            Ueta0 = Ueta;
            count = count + 1;
        end
        y = Ueta;
        yd = O.U * O.convolution_xd( O.U.' * (O.Ft - S_eta) );
        
    case 'qp'
        error('Reformulation is not implemented for quasiperiodic systems!')
end
end