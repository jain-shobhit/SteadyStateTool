classdef SSR < handle
    %SSR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        M               % Mass matrix
        K               % Stiffness matrix
        C               % Damping matrix (assuming proportional damping)
        T               % Time Period of Oscillation
        Omega           % Frequency basis vector (for quasi-periodic oscillation)
        U               % Matrix of (undamped) eigenvectors
        S               % function handle for the nonlinearity
        DS              % function handle for the Jacobian of the nonlinearity
        f               % function handle for the forcing f(t,T)
        Df              % function handle for derivative of the forcing w.r.t time period T
        order = 1       % order if integration: Newton Cotes is used (existing steps are refined)
        n_steps = 50    % number of intervals in the time domain: 50 by default
        %     end
        %
        %     properties (Access = private)
        alpha           % constants defined in (27)
        omega           % constants defined in (27)
        beta            % constants defined in (27)
        gamma           % constants defined in (27)
        omega0          % undamped natural frequencies
        o               % indices of overdamped modes
        c               % indices of critically damped modes
        u               % indices of underdamped modes
        n               % number of DOFs
        t               % time vector 0:dt:T
        dt              % node spacing in the time mesh (T/nt)
        Lvec            % Green's function for displacements (calculated over -T:dt:T)
        DLvec           % Derivative dLdT
        Jvec               % Green's function for velocities (calculated over -T:dt:T)
        ConvMtx         % convolution matrix for displacements
        isupdated       % check if different quantitities are updated according to T
        weights         % Weights vector for higher order integration
        nt              % number of nodes in the time mesh
        Ft              % forcing evaluated over time mesh
    end
    
    methods
        function O = SSR(M,C,K)
            % Basic system properties
            O.M = M;
            O.C = C;
            O.K = K;
            O.n = length(M);
            
            [V, ~] = eig(full(K),full(M));
            mu = diag(V.' * M * V);
            O.U = V * diag( 1./ sqrt(mu) ); % mass normalized VMs
            O.omega0 = sqrt(diag((O.U.' * K * O.U)));
            zeta = diag((O.U.' * C * O.U))./ (2*O.omega0);
            O.o = find(zeta>1);
            O.c = find(zeta==1);
            O.u = find(zeta<1);
            lambda2 = (-zeta - sqrt(zeta.^2 - 1)).* O.omega0;
            O.alpha = real(lambda2);
            O.omega = abs(imag(lambda2));
            O.beta = O.alpha + O.omega;
            O.gamma = O.alpha - O.omega;
        end
        
        function set.T(O,T)
            if isempty(O.T) || O.T ~= T
                O.T = T;
                update_t(O);
            end
        end
        
        function set.Omega(O,Omega)
            if isempty(O.Omega) || O.Omega ~= Omega
                O.Omega = Omega;
                not_updated(O);
            end
        end
        
        function not_updated(O)
            update.L = false;
            update.J = false;
            update.CML = false;
            update.DL = false;
            update.Ft = false;
            O.isupdated = update;
        end
        
        function set.order(O,order)
            if O.order ~= order
                O.order = order;
                update_t(O);
            end
        end
        
        function set.n_steps(O,n_steps)
            O.n_steps = n_steps;
            update_t(O);
        end
        
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
        
        function Lvec = get.Lvec(O)
            if O.isupdated.L
                Lvec = O.Lvec;
            else
                update_Lvec(O);
                Lvec = O.Lvec;
            end
        end
        
        function Ft = get.Ft(O)
            if O.isupdated.Ft
                Ft = O.Ft;
            else
                update_Ft(O);
                Ft = O.Ft;
            end
        end
        
        function update_Ft(O)
            O.Ft = O.f(O.t,O.T);
            O.isupdated.Ft = true;
        end
        
        function update_Lvec(O)
            t1 = -O.T:O.dt:O.T; % padding Green's function to ensure correct convolution
            O.Lvec = L(O,t1,O.T);
            O.isupdated.L = true;
        end
        
        function L = L(O,t,T)
            % t must be a row vector
            L = zeros(O.n,length(t));
            % underdamped
            om = repmat(O.omega(O.u), size(t));
            exp_alpha_T = repmat(exp(O.alpha(O.u) * T),size(t));
            exp_alpha_t_om = exp(O.alpha(O.u) * t)./om;
            exp_alpha_t_T_om = exp(O.alpha(O.u) * (t + T))./om;
            hu = repmat(SSR.h(t),size(O.omega(O.u)));
            denominator = repmat(1 + exp(2*O.alpha(O.u)*T) -...
                2 * exp(O.alpha(O.u)*T).* cos(O.omega(O.u)*T), size(t));
            L(O.u,:) = (exp_alpha_t_T_om.* ( sin(O.omega(O.u)*(O.T+t)) - ...
                exp_alpha_T.* sin(O.omega(O.u) * t) ) )./denominator ...
                + hu.* exp_alpha_t_om.* sin(O.omega(O.u)*t);
            
            % critically damped
            if ~isempty(O.c)
                t_c = repmat(t,size(O.c));
                exp_alpha_t = exp(O.alpha(O.c) * t);
                exp_alpha_t_T = exp(O.alpha(O.c) * (t + T));
                one_min_exp_alpha_T = 1 - repmat(exp(O.alpha(O.c) * T),size(t));
                hc = repmat(SSR.h(t),size(O.omega(O.c)));
                L(O.c,:) = exp_alpha_t_T.*(one_min_exp_alpha_T.*t_c + ...
                    T)./(one_min_exp_alpha_T.^2) + hc.*exp_alpha_t.*t_c;
            end
            
            % overdamped
            if ~isempty(O.o)
                ho = repmat(SSR.h(t),size(O.omega(O.o)));
                beta_min_gamma = repmat(O.beta(O.o) - O.beta(O.o),size(t));
                exp_beta_t_T = exp(O.beta(O.o) * (t + T));
                exp_gamma_t_T = exp(O.gamma(O.o) * (t + T));
                exp_beta_T = repmat(exp(O.beta(O.o) * T),size(t));
                exp_gamma_T = repmat(exp(O.gamma(O.o) * T),size(t));
                exp_beta_t = exp(O.beta(O.o) * t);
                exp_gamma_t = exp(O.gamma(O.o) * t);
                
                L(O.o,:) = ( (exp_beta_t_T./(1 - exp_beta_T) - ...
                    exp_gamma_t_T./(1 - exp_gamma_T)) + ho.*(exp_beta_t - ...
                    exp_gamma_t) )./ beta_min_gamma;
            end
        end
        
        function Jvec = get.Jvec(O)
            if O.isupdated.J
                Jvec = O.Jvec;
            else
                update_Jvec(O);
                Jvec = O.Jvec;
            end
        end
        
        function update_Jvec(O)
            t1 = -O.T:O.dt:O.T; % padding Green's function to ensure correct convolution
            O.Jvec = J(O,t1,O.T);
            O.isupdated.J = true;
        end
        
        function J = J(O,t,T)
            % t must be a row vector
            J = zeros(O.n,length(t));
            % underdamped
            om = repmat(O.omega(O.u), size(t));
            al = repmat(O.alpha(O.u), size(t));
            
            exp_alpha_T = repmat(exp(O.alpha(O.u) * T),size(t));
            exp_alpha_t_om = exp(O.alpha(O.u) * t)./om;
            exp_alpha_t_T_om = exp(O.alpha(O.u) * (t + T))./om;
            hu = repmat(SSR.h(t),size(O.omega(O.u)));
            denominator = repmat(1 + exp(2*O.alpha(O.u)*T) -...
                2 * exp(O.alpha(O.u)*T).* cos(O.omega(O.u)*T), size(t));
            J(O.u,:) = (exp_alpha_t_T_om.* ( om.*(cos(O.omega(O.u)*(t+T)) - ...
                exp_alpha_T.* cos(O.omega(O.u)*t)) +...
                al.*(sin(O.omega(O.u)*(O.T+t)) - ...
                exp_alpha_T.* sin(O.omega(O.u)*t))))./denominator ...
                + hu.* (exp_alpha_t_om.* (om.*cos(O.omega(O.u)*t) + ...
                al.*sin(O.omega(O.u)*t)));
            
            % critically damped
            if ~isempty(O.c)
                t_c = repmat(t,size(O.c));
                al = repmat(O.alpha(O.c), size(t));
                exp_alpha_t = exp(O.alpha(O.c) * t);
                exp_alpha_t_T = exp(O.alpha(O.c) * (t + T));
                one_min_exp_alpha_T = 1 - repmat(exp(O.alpha(O.c) * T),size(t));
                hc = repmat(SSR.h(t),size(O.omega(O.c)));
                J(O.c,:) = exp_alpha_t_T.*(one_min_exp_alpha_T.*(1 + al.*t_c) +...
                    al.*T)./(one_min_exp_alpha_T.^2) + hc.*exp_alpha_t.*(1 + ...
                    al.*t_c) ;
            end
            
            % overdamped
            if ~isempty(O.o)
                ho = repmat(SSR.h(t),size(O.omega(O.o)));
                beta_min_gamma = repmat(O.beta(O.o) - O.beta(O.o),size(t));
                bet = repmat(O.beta(O.o),size(t));
                gam = repmat(O.gamma(O.o),size(t));
                exp_beta_t_T = exp(O.beta(O.o) * (t + T));
                exp_gamma_t_T = exp(O.gamma(O.o) * (t + T));
                exp_beta_T = repmat(exp(O.beta(O.o) * T),size(t));
                exp_gamma_T = repmat(exp(O.gamma(O.o) * T),size(t));
                exp_beta_t = exp(O.beta(O.o) * t);
                exp_gamma_t = exp(O.gamma(O.o) * t);
                
                J(O.o,:) = ( (bet.*exp_beta_t_T./(1 - exp_beta_T) - ...
                    gam.*exp_gamma_t_T./(1 - exp_gamma_T)) + ...
                    ho.*(bet.*exp_beta_t - gam.*exp_gamma_t) )./ beta_min_gamma;
            end
        end
        
        function CM = get_ConvMtx(O)
            if O.isupdated.CML
                CM = O.ConvMtx;
            else
                update_ConvMtx(O);
                CM = O.ConvMtx;
            end
        end
        
        function update_ConvMtx(O)
            % This function computes the derivative of the convolution map.
            N = O.nt * O.n;
            I = zeros(O.n*O.nt^2,1);
            J = zeros(O.n*O.nt^2,1);
            val = zeros(O.n*O.nt^2,1);
            for j = 1:O.n
                ML = convmtx(O.Lvec(j,:).', O.nt);
                ML = ML(O.nt : 2*O.nt-1,:);
                I((j-1)*O.nt^2 + 1 : j*O.nt^2) = repmat(j:O.n:N,1,O.nt);
                col = repmat(j:O.n:N,O.nt,1);
                J((j-1)*O.nt^2 + 1 : j*O.nt^2) = col(:);
                val((j-1)*O.nt^2 + 1 : j*O.nt^2) = ML(:);
            end
            CM = sparse(I,J,val,N,N);
            if O.order
                w = repmat(O.weights,O.n,1);
                O.ConvMtx = CM * spdiags(w(:),0,sparse(N,N) );
            else
                O.ConvMtx = CM * O.weights;
            end
            O.isupdated.CML = true;
        end
        
        function eta = convolution_x(O,phi)
            % convolution with L to obtain velocities
            eta = zeros(size(phi));
            for j = 1:O.n
                eta(j,:) = conv(O.weights.* phi(j,:),O.Lvec(j,:),'same');
            end
        end
        
        function etad = convolution_xd(O,phi)
            % convolution with J to obtain velocities
            etad = zeros(size(phi));
            for j = 1:O.n
                etad(j,:) = conv(O.weights.* phi(j,:),O.Jvec(j,:),'same');
            end
        end
        
        function F_P = F_P(O,x)
            % function to evaluate the zero function for periodic forcing
            % eta  - modal unknowns
            % S - function handle to the nonlinearity
            % f - external forcing
            S_array = SSR.evaluate_fun_over_array(O.S,x,false);
            F_P = x + O.U * O.convolution_x(O.U.'*(S_array-O.Ft));
        end
        
        function DF_P = DF_P(O,x)
            % function to evaluate derivative of the zero function for
            % periodic forcing, currently quite expensive
            DSC = cell(O.nt,1);
            BU = cell(O.nt,1);
            for j = 1:O.nt
                DSC{j} = O.U.' * O.DS( x(:,j) ) ;
                BU{j} = O.U;
            end
            D = spblkdiag(DSC{:});
            BU = spblkdiag(BU{:});
            DF_P = speye( O.nt * O.n ) + BU * get_ConvMtx(O) * D;
        end
        
        function [x,xd] = LinearResponse(O)
            % Compute modal forcing
            phi = O.U.' * O.Ft;
            % Compute modal response by convolution with Green's function
            eta = convolution_x(O,phi);
            etad = convolution_xd(O,phi);
            x = O.U*eta;
            xd = O.U*etad;
        end
        
        function [x, xd] = Picard(O,varargin)
            % perform Picard iteration starting with optional initial guess
            % optional maximum number of iterations, and optional tolerance
            
            [x0,tol,maxiter] = parse_iteration_inputs(O,varargin);
            %% Initialization
            count = 1;
            r = 1;
            
            %% Iteration
            while r>tol
                % Evaluate the map G_P
                S_array = SSR.evaluate_fun_over_array(O.S,x0,false);
                eta = O.convolution_x( O.U.' * (O.Ft - S_array) );
                x = O.U * eta;
                % convergence check
                r = O.dt * norm(x-x0,Inf)/norm(x0,Inf);
                disp(['Iteration ' num2str(count) ': ' 'relative change = ' num2str(r)])
                if count>maxiter || isnan(r)
                    error('Picard iteration did not converge')
                end
                x0 = x;
                count = count + 1;
            end
            xd = O.U * O.convolution_xd( O.U.' * (O.Ft - S_array) );
        end
        
        function [x, xd] = NewtonRaphson(O,varargin)
            % perform Newton--Raphson iteration starting with an optional
            % initial guess, optional maximum number of iterations and
            % optional tolerance
            [x0,tol,maxiter] = parse_iteration_inputs(O,varargin);
            
            count = 1;
            r=1;
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
        end
        
        function [x0,tol,maxiter] = parse_iteration_inputs(O,inputs)
            % this method used for parsing interation inputs
            defaultTol = 1e-6;
            defaultInit = O.LinearResponse();
            defaultmaxiter = 20;
            p = inputParser;
            addOptional(p,'tol',defaultTol)
            addOptional(p,'init',defaultInit)
            addOptional(p,'maxiter',defaultmaxiter)
            parse(p,inputs{:});
            x0 = p.Results.init;
            maxiter = p.Results.maxiter;
            tol = p.Results.tol;
        end
        
        function [T_array, sol] = sequential_continuation(O,T_range,varargin)
            % this method performs sequential continuation of periodic
            % response assuming the time-period (T) of the solution to be
            % the continuation parameter. Warning: this scheme cannot
            % capture fold points in the frequency response curves.
            
            T0 = T_range(1); % the time period for the initial solution
            O.T = T0;
            x0 = LinearResponse(O);
            n_cont_steps = 100; % number of continuation steps
            
            dT = diff(T_range)/n_cont_steps;
            T_array = T_range(1):dT:T_range(2);
            sol = cell(1,n_cont_steps);
            try_Picard = true;
            
            for j = 1:n_cont_steps
                T0 = T_array(j);
                O.T = T0;
                
                if isempty(x0)
                    x0 = LinearResponse(O);
                end
                
                if try_Picard
                    try
                        [x0, xd0] = Picard(O,'init',x0);
                    catch
                        disp('Switching to Newton--Raphson iteration')
                        try_Picard = false;
                        [x0, xd0] = NewtonRaphson(O,'init',x0);
                    end
                else
                    [x0, xd0] = NewtonRaphson(O,'init',x0);
                end                                
                sol{j} = [x0;xd0];                
            end
            
        end
    end
    
    methods(Static)
        function w = NewtonCotes(order)
            % returns the weights corresponding to the appropriate order for performing
            % Newton-Cotes numerical integration on a uniformly spaced grid
            switch order
                case 1
                    w = [1 1]/2;
                case 2
                    w = [1 4 1]/3;
                case 3
                    w = [1 3 3 1]*(3/8);
                case 4
                    w = [7 32 12 32 7]*(2/45);
            end
        end
        
        function h = h(t)
            h = t>0; % use overloading of the logical operator
        end
        
        function F = evaluate_fun_over_array(f,x,vectorized)
            % accepts single argument functions and evaluates them over
            % array x
            if vectorized
                F = f(x);
            else
                F0 = f(x(:,1));
                F = zeros(size(F0,1),size(x,2));
                F(:,1) = F0;
                for j = 2:size(x,2)
                    F(:,j) = f(x(:,j));
                end
            end
        end
    end
end

