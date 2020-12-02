function update_Gvec(O)
% Computing Green's function over grid [0,T]
G1 = G(O,O.t,O.T);
% Padding Green's function, i.e., using values over grid [-T,T] 
% to ensure correct convolution (using the fact that Green's function 
% is T periodic) 
O.Gvec = [G1(:,1:end-1) G1];
O.isupdated.G = true;
end