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
