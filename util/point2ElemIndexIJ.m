function index = point2ElemIndexIJ(x,mesh)
%
% given a particle position x, and a background mesh 'mesh'
% find the index (i,j) of the element/cell contains 'x'.
%
% Inspired by Uintah-MPM implementation.
%  
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% February 2014.
% June 2014.


deltax = mesh.deltax;
deltay = mesh.deltay;

xmin   = min(mesh.node(:,1));
ymin   = min(mesh.node(:,2));

index.i = floor((x(1)-xmin)/deltax) + 1;
index.j = floor((x(2)-ymin)/deltay) + 1;

% handle objects cut by the boundaries of the grid

if index.i > mesh.numx
  index.i = mesh.numx;
end

if index.i <= 0
  index.i = 1;
end

if index.j > mesh.numy
  index.j = mesh.numy;
end

if index.j <= 0
  index.j = 1;
end


