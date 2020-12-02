function DF_P = DF_P_reformulated(O,UTz)
% function to evaluate derivative of the zero function for
% periodic forcing, reformulated setting
DSC = cell(O.nt,1);
for j = 1:O.nt
    DSC{j} = sparse(O.U.' * O.DS( UTz(:,j) ) * O.U );
end
D = blkdiag(DSC{:});
A = get_AA(O);

p = gcp('nocreate');        % check if there is a parallel pool

% only distribute arrays for large enough systems
if (O.n > 150) && (~isempty(p)) 
    A = distributed(A);
    D = distributed(D);
    spmd
        DD = D * A;
    end
    DF_P = distributed.speye(size(DD)) + DD;
    DF_P = full(DF_P);      % DF_P is always full
else
    DD = D * A;
    DF_P = speye(size(DD)) + DD;
    DF_P = full(DF_P);      % DF_P is always full
end

end