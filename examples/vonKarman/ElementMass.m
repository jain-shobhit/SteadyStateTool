function M = ElementMass(prop)

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

