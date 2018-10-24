function update_Lvec(O)
t1 = -O.T:O.dt:O.T; % padding Green's function to ensure correct convolution
O.Lvec = L(O,t1,O.T);
O.isupdated.L = true;
end