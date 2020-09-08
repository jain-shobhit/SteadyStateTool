function DF_P = DF_P(O,x)
% function to evaluate derivative of the zero function for
% periodic forcing, currently quite expensive
DSC = cell(O.nt,1);
BU = cell(O.nt,1);
for j = 1:O.nt
    DSC{j} = O.U.' * O.DS( x(:,j) ) ;
    BU{j} = O.U;
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
end