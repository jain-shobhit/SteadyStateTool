function DF_P = DF_P(O,x)
% function to evaluate derivative of the zero function for
% periodic forcing, currently quite expensive
DSC = cell(O.nt,1);
switch O.sys_order
    case 'first'
        BU = cell(O.nt,1);
        z = x;
        for j = 1:O.nt
            DSC{j} = O.W' * O.DR( z(:,j) ) ;
            BU{j} = O.V;
        end
        if verLessThan('matlab','9.7.0')
            D = spblkdiag(DSC{:});
            BU = spblkdiag(BU{:});
        else
            D = blkdiag(DSC{:});
            BU = blkdiag(BU{:});
        end
        DD = BU * get_ConvMtx(O) * D;
        DF_P = speye( size(DD) ) + DD;
    case 'second'
        B = size(O.U);
        for j = 1:O.nt
            if B(1) == B(2)
                DSC{j} = sparse(O.U.' * O.DS( x(:,j) ) );
            else
                DSC{j} = sparse(O.U.' * O.DS( O.U * x(:,j) ) * O.U);
            end
        end
        D = blkdiag(DSC{:});
        C = get_ConvMtx(O);         % C is already full here

        p = gcp('nocreate');        % check if there is a parallel pool

        % only distribute arrays for large enough systems
        if (O.n > 150) && (~isempty(p)) 
            C = distributed(C);
            D = distributed(D);
            spmd
                DD = C * D;
            end
            DF_P = distributed.speye( size(DD) ) + DD;
        else
            DD = C * D;
            DF_P = speye(size(DD)) + DD;
        end
end
end