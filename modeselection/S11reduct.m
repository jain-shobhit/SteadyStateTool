function S = S11reduct(SS,X)
% This function converts the quadratic nonlinearity coefficients from
% their general form to the desired form in modal coordinates
    U = SS.U;
    U2 = SS.U2;
    B = size(U2);
    S = cell(1,B(2));
    n = B(2)+SS.n;
    
    for j = 1:B(2)
        S{j} = zeros(SS.n,SS.n);
        for i = 1:n
            S{j} = S{j} + U2(i,j) * U.' * X{i} * U;
        end
    end    
end