function DS = nonsmooth_nonlinearity_derv(x,alpha,beta)
modx = abs(x);    
if modx>beta
    DS = alpha;
else
    DS = 0;
end

end