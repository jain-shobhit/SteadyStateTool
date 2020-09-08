function M = ElementMass(prop)
% The finite element model is taken from the following article: 
% Jain, S., Tiso, P., & Haller, G. (2018). Exact nonlinear model reduction 
% for a von Kármán beam: slow-fast decomposition and spectral submanifolds. 
% Journal of Sound and Vibration, 423, 195–211. 
% https://doi.org/10.1016/J.JSV.2018.01.049

rho = prop.rho;
A   = prop.A;
l   = prop.dx;


M = sparse(6,6);

M([1 4],[1 4]) = ((rho*A*l)/6)*[2 1;
                                        1 2];
                              
M([2 3 5 6],[2 3 5 6]) = ((rho*A*l)/420)*[156      22*l       54      -13*l;
                                      22*l     4*l^2      13*l    -3*l^2;
                                      54       13*l       156     -22*l;  
                                      -13*l     -3*l^2     -22*l   4*l^2];

