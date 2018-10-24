function S = nonsmooth_nonlinearity(x,alpha,beta)
modx = abs(x);    
if modx>beta
    S = alpha*sign(x) * (modx - beta);
else
    S = 0;
end

end