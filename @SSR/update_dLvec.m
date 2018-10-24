function update_dLvec(O)
t1 = -O.T:O.dt:O.T; % padding Green's function to ensure correct convolution
O.DLvec = dLdT(O,t1,O.T);
O.isupdated.dL = true;
end