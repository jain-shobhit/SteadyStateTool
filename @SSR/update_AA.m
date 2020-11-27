function update_AA(O)
% constructs the matrix equivalent of the A operation in DF
Lint = get_Lint(O);
N = O.nt * O.n;
I = zeros(O.n*O.nt^2,1);
J = zeros(O.n*O.nt^2,1);
val = zeros(O.n*O.nt^2,1);
for j = 1:O.n
    ML = squeeze(Lint(j,:,:)).';
    I((j-1)*O.nt^2 + 1 : j*O.nt^2) = repmat(j:O.n:N,1,O.nt);
    col = repmat(j:O.n:N,O.nt,1);
    J((j-1)*O.nt^2 + 1 : j*O.nt^2) = col(:);
    val((j-1)*O.nt^2 + 1 : j*O.nt^2) = ML(:);
end
A = sparse(I,J,val,N,N);
O.AA = A;
O.isupdated.AA = true;
end