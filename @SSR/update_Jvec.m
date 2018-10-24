function update_Jvec(O)
t1 = -O.T:O.dt:O.T; % padding Green's function to ensure correct convolution
O.Jvec = J(O,t1,O.T);
O.isupdated.J = true;
end