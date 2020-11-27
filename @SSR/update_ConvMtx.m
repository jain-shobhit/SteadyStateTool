function update_ConvMtx(O)
% This function computes the derivative of the convolution map.

switch O.sys_order
    case 'first'
        ntN = O.nt * O.N;
        I = zeros(O.N*O.nt^2,1);
        J = zeros(O.N*O.nt^2,1);
        val = zeros(O.N*O.nt^2,1);
        for j = 1:O.N
            MG = convmtx(O.Gvec(j,:).', O.nt);
            MG = MG(O.nt : 2*O.nt-1,:);
            I((j-1)*O.nt^2 + 1 : j*O.nt^2) = repmat(j:O.N:ntN,1,O.nt);
            col = repmat(j:O.N:ntN,O.nt,1);
            J((j-1)*O.nt^2 + 1 : j*O.nt^2) = col(:);
            val((j-1)*O.nt^2 + 1 : j*O.nt^2) = MG(:);
        end
        CM = sparse(I,J,val,ntN,ntN);
        if O.order
            w = repmat(O.weights,O.N,1);
            O.ConvMtx = CM * spdiags(w(:),0,sparse(ntN,ntN) );
        else
            O.ConvMtx = CM * O.weights;
        end
        O.isupdated.CML = true;
        
    case 'second'
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
        
        % This is here so that BU*CM happens once per cont. step rather
        % than once per iteration step
        B = size(O.U);
        if B(1) == B(2)
            BU = cell(O.nt,1);
            for j = 1:O.nt
                BU{j} = sparse(O.U); 
            end
            BU = blkdiag(BU{:});
            CM = BU * CM;
        end
        
        CM = full(CM); % at this point CM is full

        if O.order
            w = repmat(O.weights,O.n,1);
            O.ConvMtx = CM * spdiags(w(:),0,sparse(N,N) );
        else
            O.ConvMtx = CM * O.weights;
        end
        O.isupdated.CML = true;
end
