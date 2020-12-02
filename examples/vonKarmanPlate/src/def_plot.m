function []=def_plot(u0,scale,nodes,elements)
%
% []=def_plot(u0,scale,nodes,elements)

u0=reshape(u0,6,size(nodes,1))';

% for i=1:size(u0,1)
%     magn(i)=(norm(u0(i,1:3)))^2;
% end
%
% maxdispl=max(magn);

xmax = max(nodes(:,2)); xmin = min(nodes(:,2));
ymax = max(nodes(:,3)); ymin = min(nodes(:,3));
zmax = max(nodes(:,4)); zmin = min(nodes(:,4));

d = sqrt((xmax-xmin)^2+(ymax-ymin)^2+(zmax-zmin)^2);

u=u0(:,1);
v=u0(:,2);
w=u0(:,3);
displ=sqrt(u.^2+v.^2+w.^2);
maxdispl=max(displ);

sf=scale*d/maxdispl;
if scale == 'TRUE'
    sf = 1;
end

% figure
trisurf(elements(:,2:4),nodes(:,2)+sf*u,nodes(:,3)+sf*v,nodes(:,4)+sf*w,0*displ)
hold on
%trimesh(elements(:,2:4),nodes(:,2),nodes(:,3),nodes(:,4),0*displ)
axis equal
%     axis off
xlabel('x (mm)')
ylabel('y (mm)')
% hold off
%shading interp
% scale=0.5*(100)/maxdispl;
%
% figure
% trisurf(elements(:,2:end),nodes(:,2)+scale*u0(:,1),nodes(:,3)+scale*u0(:,2),nodes(:,4)+scale*u0(:,3))

% patch('Vertices',[nodes(:,2)+scale*u0(:,1) nodes(:,3)+scale*u0(:,2) nodes(:,4)+scale*u0(:,3)],'Faces',elements(:,2:end),...
%     'FaceColor','interp')
%axis([-100 100 -100 600 -50 150])
% view(-90,0)
% for the curved panel
% axis([0.1676   50.0000         0   50.0000   50   100]);
% for the short cantilever
% axis([-29   51   -6  180         0  300.0])
% for the hole plate
% axis([-80.3888   40.3921  -30.0000   30.0000   -5.5659    6.5225])

% axis([-10   470     -10    30     -20    30])
% hold off
