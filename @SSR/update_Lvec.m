function update_Lvec(O)
% Computing Green's function over grid [0,T]
L1 = L(O,O.t,O.T);
% Padding Green's function, i.e., using values over grid [-T,T] 
% to ensure correct convolution (using the fact that Green's function 
% is T periodic) 
O.Lvec = [L1(:,1:end-1) L1];
O.isupdated.L = true;
end