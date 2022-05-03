function e = point2ElemIndex3D(x,mesh)
%
% given a particle position x, and a background mesh 'mesh'
% find the index of the element/cell contains 'x'.
%
% Vinh Phu Nguyen
% Monash University, Victoria, Australia
% April 2016.


deltax = mesh.deltax;
deltay = mesh.deltay;
deltaz = mesh.deltaz;

xmin   = min(mesh.node(:,1));
ymin   = min(mesh.node(:,2));
zmin   = min(mesh.node(:,3));

e = floor((x(1)-xmin)/deltax) + 1 + mesh.numx*floor((x(2)-ymin)/deltay)+...
  + mesh.numx*mesh.numy*floor((x(3)-zmin)/deltaz);
