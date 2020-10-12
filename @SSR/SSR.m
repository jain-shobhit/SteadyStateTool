classdef SSR < handle
    %SSR (SteadyStateResponse) class implements the integral equations
    %approach to the computation of the steady-state response in
    %(quasi)periodically-forced dynamical systems.
    
    % The theory of computation (without model reduction) is decribed in
    % the following article:
    %
    
    % The theory for model reduction is described in the following article:
    % G. Buza, S. Jain, G. Haller, Using Spectral Submanifolds for Optimal
    % Mode Selection in Model Reduction, (2020) Preprint available on arXiv.org
    
    properties
        
        sys_order       % order in which system is set up, can be 'first', 'second'
        
       
        %First order
        A               % linear System matrix Bdot(z) = Az +F(t)
        B               % System matrix 
        V               % Matrix of eigenvectors
        W
        R               % function handle for the nonlinearity, takes the 2n dim vector z as input
        DR              % function handle for the Jacobian of the nonlinearity
        F               % function handle for the forcing F(t,T)
        DF              % function handle for derivative of the forcing w.r.t time period T
        
        %Second order
        M               % Mass matrix
        K               % Stiffness matrix
        C               % Damping matrix (assuming proportional damping)
        U               % Matrix of (undamped) eigenvectors
        S               % function handle for the nonlinearity
        DS              % function handle for the Jacobian of the nonlinearity
        f               % function handle for the forcing f(t,T)
        Df              % function handle for derivative of the forcing w.r.t time period T

        % Common parameters
        T               % Time Period of Oscillation
        Omega           % Frequency basis vector (for quasi-periodic oscillation)
        order = 1       % order if integration: Newton Cotes is used (existing steps are refined)
        n_steps = 50    % number of intervals in the time domain: 50 by default
        type = 'p'      % 'p': calculations performed for periodic solution (time domain)
                        % 'qp' : calculations performed for periodic solution (time domain)
        %     end
        
        
        %     properties (Access = private)
        alpha           % constants defined in (27)
        omega           % constants defined in (27) 
                        % - eigenfrequencies of linear system in case of
                        %   first order implementation
        beta            % constants defined in (27)
        gamma           % constants defined in (27)
        omega0          % undamped natural frequencies
        o               % indices of overdamped modes
        c               % indices of critically damped modes
        u               % indices of underdamped modes
        n               % number of DOFs
        t               % time vector 0:dt:T
        dt              % node spacing in the time mesh (T/nt)
        ConvMtx         % convolution matrix for displacements
        isupdated       % check if different quantitities are updated according to T
        weights         % Weights vector for higher order integration
        nt              % number of nodes in the time mesh
        Ft              % forcing evaluated over time mesh
        nh              % number of harmonics in of the base frequency Omega
        kappa_set       % the set of frequency-indices used during Fourier transformation
        theta_set       % the set of torus coordinates over which the solution is approximated used during Fourier transformation
        
        % First order implementation
        Gvec            % Green's function first order (calculated over -T:dt:T)
        Hmat            % the global H matrix containing the Green's function (amplification factors) for freq domain analysis
        F_theta         % function handle for forcing in torus coordinates
        F_kappa         % Fourier coefficients for the external forcing
        N                % Phase space size        
        
        % Second order implementation
        Lvec            % Green's function for displacements (calculated over -T:dt:T)
        DLvec           % Derivative dLdT
        Lint            % precomputed integrals of Green's functions
        Jvec            % Green's function for velocities (calculated over -T:dt:T)
        Qmat            % the global Q matrix containing the Green's function (amplification factors) for freq domain analysis
        f_theta         % function handle for forcing in torus coordinates
        f_kappa         % Fourier coefficients for the external forcing
        lambda          % contains eigenvalues of the system
        
        m               % the discretization along each dimension of the Torus
        E               % The Fourier transform exponential matrix
        Einv            % The inverse Fourier transform exponential matrix
        Evinv           % The inverse Fourier transform exponential matrix
        Ekron           % The Fourier transform exponential matrix kronecker product with Identity
        Einvkron        % The inverseFourier transform exponential matrix kronecker product with Identity
        mode_choice     % Selection of modes specified by the user
        
        % required only for the SSM computation in the mode selection
        % procedure
        zeta            % The damping coefficients of the master modes
        omega2          % The eigenfrequencies of the enslaved modes
        zeta2           % The damping coefficients of the enslaved modes
        K1              % Stiffness matrix of the enslaved modes (in modal coordinates)
        C1              % Damping matrix of the enslaved modes (in modal coordinates)
        U2              % Enslaved eigenmodes
        % storage for some of these
        U2store
        zeta2store
        omega2store
        SSMdone         % Boolean that ensures that the eigenvalues are only computed once
    end
    
    methods
        % Initiate an element of the SSR class by calling this function with
        % a struct 'System' as input.
        % System contains system matrices and string stating if system is
        % set up in first or second order form
        function O = SSR(System,varargin)
            
            switch System.sys_order
                case 'first'
                    % Basic system properties
                    O.sys_order = 'first';
                    O.A = System.A;
                    O.B = System.B;
                    
                    O.n = length(O.A)/2;
                    O.N = length(O.A);
                    [V, lambda, W] = eig(full(O.A),full(O.B));
                    mu = diag(W' * O.B * V);
                    O.V = V * diag( 1./ sqrt(mu) ); % mass normalized VMs
                    O.W = (W*diag(1./(sqrt(mu)')));%O.V.';
                    
                    O.lambda = diag(lambda); % is now a vector
                    omega  = abs(imag(O.lambda)); % natural frequencies
                    O.omega = unique(omega);
                    
                case 'second'
                    % Basic system properties
                    O.sys_order = 'second';
                    O.M = System.M;
                    O.C = System.C;
                    O.K = System.K;
                    
                    if nargin == 1
                        n = length(O.M);
                        O.mode_choice = 1:n;
                        [VV, dd] = eig(full(O.K),full(O.M));
                    else
                        O.mode_choice = varargin{1};
                        [VV, dd] = eigs(O.K,O.M,max(O.mode_choice),'SM');
                    end
                    O.n = length(O.mode_choice);
                    
                    dd = diag(dd);
                    [~, ind] = sort(dd);
                    VV = VV(:,ind);
                    V = VV(:,O.mode_choice);
                    mu = diag(V.' * O.M * V);
                    O.U = V * diag( 1./ sqrt(mu) ); % mass normalized VMs
                    O.omega0 = sqrt(diag((O.U.' * O.K * O.U)));
                    O.K1 = O.U.' * O.K * O.U;
                    O.zeta = diag((O.U.' * O.C * O.U))./ (2*O.omega0);
                    O.C1 = O.U.' * O.C * O.U;
                    O.o = find(O.zeta>1);
                    O.c = find(O.zeta==1);
                    O.u = find(O.zeta<1);
                    lambda1 = (-O.zeta + sqrt(O.zeta.^2 - 1)).* O.omega0;
                    lambda2 = (-O.zeta - sqrt(O.zeta.^2 - 1)).* O.omega0;
                    O.alpha = real(lambda2);
                    O.omega = abs(imag(lambda2));
                    O.beta = lambda1;
                    O.gamma = lambda2;
                    O.SSMdone = 0;
            end
        end  
        
        function U2 = get.U2(O)
            if ~O.SSMdone
                SSMcomps(O);
            end
            U2 = O.U2store;
        end
        
        function omega2 = get.omega2(O)
            if ~O.SSMdone
                SSMcomps(O);
            end
            omega2 = O.omega2store;
        end
        
        function zeta2 = get.zeta2(O)
            if ~O.SSMdone
                SSMcomps(O);   
            end
            zeta2 = O.zeta2store;
        end
        
        function set.T(O,T)
            if ~isequal(O.T,T)
                O.T = T;
                update_t(O);
            end
        end
        
        function set.Omega(O,Omega)
            if ~isrow(Omega)
                Omega = Omega.';
            end
            if ~isequal(O.Omega, Omega)
                O.Omega = Omega;
                not_updated(O);
            end
        end
        
        notupdated(O)
        
        function set.order(O,order)
            if O.order ~= order
                O.order = order;
                update_t(O);
            end
        end
        
        update_t(O)
        
        function set.n_steps(O,n_steps)
            O.n_steps = n_steps;
            update_t(O);
        end        
       
        function set.nh(O,nh)
            if ~isequal(O.nh, nh)
                O.nh = nh;
                not_updated(O);
            end
        end
        
        function set.m(O,m)
            if ~isequal(O.m, m)
                O.m = m;
                not_updated(O);
            end
        end
        
        function theta = get.theta_set(O)
            if ~O.isupdated.theta
                update_theta_set(O);
            end
            theta = O.theta_set;
        end
        
        update_theta_set(O)

        function kappa = get.kappa_set(O)
            if ~O.isupdated.kappa
                update_kappa_set(O);
            end
            kappa = O.kappa_set;
        end
        
        update_kappa_set(O)
        
        function Lvec = get.Lvec(O)
            if ~O.isupdated.L
                update_Lvec(O);
            end
            Lvec = O.Lvec;
        end

        function Gvec = get.Gvec(O)
            % updates greens function 
            if ~O.isupdated.G
                update_Gvec(O);
            end
            Gvec = O.Gvec;
        end        
        
        function Ft = get.Ft(O)
            if ~O.isupdated.Ft
                update_Ft(O);
            end
            Ft = O.Ft;
        end
        
        update_Ft(O)
        
        function E = get.E(O)
            if ~O.isupdated.E
                update_E(O);
            end
            E = O.E;
        end
        
        update_E(O)
        
        function f_kappa = get.f_kappa(O)
            if ~O.isupdated.f_kappa
                update_f_kappa(O);
            end
            f_kappa = O.f_kappa;
        end
        
        function F_kappa = get.F_kappa(O)
            if ~O.isupdated.F_kappa
                update_f_kappa(O);
            end
            F_kappa = O.F_kappa;
        end
        
        update_f_kappa(O)
        
        update_Lvec(O)
        
        update_Gvec(O)
                
        L = L(O,t,T)
        
        function Jvec = get.Jvec(O)
            if O.isupdated.J
                Jvec = O.Jvec;
            else
                update_Jvec(O);
                Jvec = O.Jvec;
            end
        end
               
        update_Jvec(O)
        
        J = J(O,t,T)
        
        CM = get_ConvMtx(O)
        
        update_ConvMtx(O)
        
        update_dLvec(O)
        
        L = dLdT(O,t,T)
        
        function Q = get.Qmat(O)
            if ~O.isupdated.Q
                update_Qmat(O);
            end
            Q = O.Qmat;
        end
        
        update_Qmat(O)
        
        function H = get.Hmat(O)
            if ~O.isupdated.H
                update_Hmat(O);
            end
            H = O.Hmat;
        end
        
        update_Hmat(O)
        
        eta = convolution_x(O,phi)
        
        etad = convolution_xd(O,phi)
        
        eta = convolution_dT(O,phi)
        
        F_P = F_P(O,x)
        
        F_Q = F_Q(O,x)
        
        DF_P = DF_P(O,x)
        
        DF_Q = DF_Q(O,x)
        
        [x, xd] = LinearResponse(O)
        
        [x, xd] = Picard(O,varargin)
        
        [x, xd] = NewtonRaphson(O,varargin)
        
        [x0,tol,maxiter] = parse_iteration_inputs(O,inputs)
        
        [OMEGA, SOL, PICARD] = sequential_continuation(O,range,varargin)
               
        W = SSM2(O,S)
        
        [U2,omega2,zeta2] = SSMcomps(O)
        
        [om, N] = FRC(O,T_range,ncontsteps)
    end
    
    methods(Static)
        w = NewtonCotes(order)
        
        h = h(t)
        
        F = evaluate_fun_over_array(f,x,vectorized)
    end
end

