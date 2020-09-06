function PlotDeformation(model,displacement,scale)

l = model.nElements*model.dx;
u = displacement(1:3:end,:); w = displacement(2:3:end,:);
dis = sqrt(u.^2 + w.^2);
maxdispl = max(dis);

sf = scale*l./maxdispl;
sf = diag(sf);
if strcmp(scale, 'TRUE')
    sf = 1;
end
% sf = scale;

x = 0:model.dx:l;
% y = zeros(size(x));
xd = x' + u*sf;
yd = w*sf;

% plot(x,y,'-k','LineWidth',2)
% hold on
plot(xd,yd)
% drawnow
% % axis equal
% plot(xd(model.nElements/4),yd(model.nElements/04), '*', 'Markersize' , 6)
