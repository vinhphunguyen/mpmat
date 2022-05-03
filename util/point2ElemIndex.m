function e = point2ElemIndex(x,mesh)
%
% given a particle position x, and a background mesh 'mesh'
% find the index of the element/cell contains 'x'.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% February 2014.
% June 2014.


deltax = mesh.deltax;
deltay = mesh.deltay;

xmin   = min(mesh.node(:,1));
ymin   = min(mesh.node(:,2));

e = floor((x(1)-xmin)/deltax) + 1 + mesh.numx*floor((x(2)-ymin)/deltay);
