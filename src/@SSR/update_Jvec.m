function update_Jvec(O)
% Computing Green's function over grid [0,T]
J1 = J(O,O.t,O.T);
% Padding Green's function, i.e., using values over grid [-T,T] 
% to ensure correct convolution (using the fact that the Green's function 
% is T periodic) 
O.Jvec = [J1(:,1:end-1) J1];
O.isupdated.J = true;
end