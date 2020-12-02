function [K2global,K3global] = RotateK2K3(K2local,K3local, T)
% this function rotates the local third order element stiffness tensor in
% to global coordinates. It is reponsible for quadratic internal forces.

nndof = length(T);   % number of dofs per element
K2 = tensor(K2local,[nndof nndof nndof]);
K3 = tensor(K3local,[nndof nndof nndof nndof]);
Q  = tensor(full(T),[nndof nndof]);
K2global = ttt(ttt(ttt(Q,K2,1,1),Q,2,1),Q,2,1);
K3global = ttt(ttt(ttt(ttt(Q,K3,1,1),Q,2,1),Q,2,1),Q,2,1);

