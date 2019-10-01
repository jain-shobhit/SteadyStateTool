function DF_Q = DF_Q(O,x)
% function to evaluate derivative of the zero function for
% quasiperiodic forcing, currently quite expensive
DSC = cell(O.nt,1);
N = O.n * size(O.kappa_set,2);
for j = 1:size(x,2)
    DSC{j} = O.DS( x(:,j) ) ;
end
if verLessThan('matlab','9.7.0')
    D = spblkdiag(DSC{:});
else
    D = blkdiag(DSC{:});
end 

DF_Q = speye(N,N) + O.Qmat * O.Ekron * D * O.Einvkron;
end