function [Geo] = Geometry(eps,a)
%% Geometry for a Von Karman Beam element aligned along the x axis (all units SI)
Geo.eps = eps;
Geo.l = 1 ; % length 
Geo.h = eps*Geo.l; % height
Geo.b = Geo.l/10; % width
Geo.A = Geo.b*Geo.h; % cross section area
Geo.I =  Geo.b*Geo.h^3/12; % second moment of cross section
Geo.a = a; % height of midpoint relative to the ends (measure of curvature of the beam)
