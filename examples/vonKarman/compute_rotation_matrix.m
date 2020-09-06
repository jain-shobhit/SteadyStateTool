function T = compute_rotation_matrix(nodes)
T = speye(6,6);

node1=nodes(1,:);
node2=nodes(2,:);
l=norm(node2-node1);
e_x = (node2-node1)/l;
R = [e_x(1) e_x(2);
    -e_x(2) e_x(1)];

T(1:2,1:2) = R;
T(4:5,4:5) = R;
end