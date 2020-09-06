function [w,VN] = sortVMsVK(V,M,D)

w = (diag(D));
[w,pos] = sort(w);
VN = V(:,pos);


mu = diag(VN'*M*VN);
VN = VN*diag(1./sqrt(mu));

end
