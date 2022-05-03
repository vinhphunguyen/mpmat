function e = point2ElemIndex1D(x,mesh)
%
% given a particle position x, and a background mesh 'mesh'
% find the index of the element/cell contains 'x'.
% For 1D problems.
%
% Vinh Phu Nguyen
% Monash University
% July 2016.
% June 2016.


deltax = mesh.deltax;
xmin   = min(mesh.node);

e = floor((x-xmin)/deltax) + 1;

if ( e > mesh.elemCount ), e = mesh.elemCount; end
