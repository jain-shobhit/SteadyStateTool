function DF_P = DF_P(O,x)
% function to evaluate derivative of the zero function for
% periodic forcing, currently quite expensive
DSC = cell(O.nt,1);
BU = cell(O.nt,1);
for j = 1:O.nt
    DSC{j} = O.U.' * O.DS( x(:,j) ) ;
    BU{j} = O.U;
end
D = spblkdiag(DSC{:});
BU = spblkdiag(BU{:});
DF_P = speye( O.nt * O.n ) + BU * get_ConvMtx(O) * D;
end