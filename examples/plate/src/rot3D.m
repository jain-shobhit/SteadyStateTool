function T=rot3D(nodes)

node1=nodes(1,:);
node2=nodes(2,:);
node3=nodes(3,:);

l=norm(node2-node1);

e_x = (node2-node1)/l;
d2  = node2-node1;
d3  = node3-node1;
e_z = cross(d2,d3)/norm(cross(d2,d3));
e_y = cross(e_z,e_x);

R = [ e_x(1) e_x(2) e_x(3);
      e_y(1) e_y(2) e_y(3);
      e_z(1) e_z(2) e_z(3)];
   
T = sparse(18,18);

T(1:3,1:3) = R;
T(4:6,4:6) = R;
T(7:9,7:9) = R;
T(10:12,10:12) = R;
T(13:15,13:15) = R;
T(16:18,16:18) = R;
