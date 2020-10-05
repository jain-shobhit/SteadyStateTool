function DF_Q = DF_Q(O,x)
% function to evaluate derivative of the zero function for
% quasiperiodic forcing, currently quite expensive
DSC = cell(O.nt,1);

switch O.sys_order
    case 'first'
        N_k = O.N * size(O.kappa_set,2);
        for j = 1:size(x,2)
            DRC{j} = O.DR( x(:,j) ) ;
        end
        if verLessThan('matlab','9.7.0')
            D = spblkdiag(DRC{:});
        else
            D = blkdiag(DRC{:});
        end
        
        DF_Q = speye(N_k,N_k) + O.Hmat * O.Ekron * D * O.Einvkron;
    case 'second'
        
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