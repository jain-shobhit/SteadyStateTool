function W = SSM2(O,S)
% A faster computation of the leading order term of the SSM

n = length(O.M);
m = O.n;
W = cell(1,n-m);

dim = (2*m)^2;

A = [zeros(m),eye(m);-O.K1,-O.C1];  % linear coefficient of master system
B0 = zeros(dim);

delta = eye(2*m);                   % kronecker delta 

% assembly of the 2D equivalent to tensor B as a part of vectorizing the
% problem (first half that is independent of enslaved mode)
for i = 1:dim
    for j = 1:dim
        q = floor((j-1)/2/m)+1;
        t = floor((i-1)/2/m)+1;
        r = j - (q-1) * 2 * m;
        s = i - (t-1) * 2 * m;
        B0(i,j) = 2 * A(r,s) * A(q,t) +...
            A(q,:) * A(:,t) * delta(r,s) +...
            A(r,:) * A(:,s) * delta(q,t);   % eq. (2.14) in thesis 
    end
end


% Computing SSM coefficient Wk for each enslaved mode
for k = 1:n-m
    B1 = zeros(dim);
    Sk = zeros(2*m);
    Sk(1:m,1:m) = S{k};     % this is R_k, nonlinear quadratic coeffs in original eq.
    
    % Computing second half of B
    for i = 1:dim
        for j = 1:dim
            q = floor((j-1)/2/m)+1;
            t = floor((i-1)/2/m)+1;
            r = j - (q-1) * 2 * m;
            s = i - (t-1) * 2 * m;
            B1(i,j) = 2 * (O.zeta2(k)) * (O.omega2(k)) *...
                (A(q,t) * delta(r,s) + A(r,s) * delta(q,t));
        end
    end    
    
    B = B0 + B1 + (O.omega2(k))^2 * eye(dim);
    Wk = -B\Sk(:);                  % quadratic SSM coeffs in vector form
    W{k} = reshape(Wk,size(Sk));    % converting SSM coeffs back to matrix form
end

end