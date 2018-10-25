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