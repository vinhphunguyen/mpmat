function [res] = generateMPForQuad(geo,ppc,mesh)
%
% Generate material points for a rectangle 'geo' where PPC (particle per cell)
% is given by 'ppc=[ppc_x ppc_y]' and background mesh is 'mesh'.
% geo is a matrix 4x2.
%
% Use built-in function 'inpolygon' to check if a point is within a
% polygon.
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% August 2015.

% bounding box of the quad
mi = min(geo);
ma = max(geo);

xmin = mi(1); ymin = mi(2);
xmax = ma(1); ymax = ma(2);

% indices of cells containing xmin, xmax, ymin, ymax
minIndex = point2ElemIndexIJ([xmin ymin],mesh);
maxIndex = point2ElemIndexIJ([xmax ymax],mesh);

iMin     = minIndex.i;
jMin     = minIndex.j;
iMax     = maxIndex.i;
jMax     = maxIndex.j;

noElems = (iMax-iMin+1) * (jMax-jMin+1);

% res.elems contains all elements cut the cirlce
res.elems = zeros(noElems,1);

dx  = mesh.deltax/(ppc(1));
dy  = mesh.deltay/(ppc(2));

volume = [];
coord  = [];
ppc(2)
ii = 1;
for i = iMin:iMax
  for j = jMin:jMax
    id  = i + mesh.numx * ( j - 1);
    res.elems(ii) = id;
    ii  = ii + 1;
    sctr = mesh.element(id,:);          %  element scatter vector
    pts  = mesh.node(sctr,:);
    x1   = pts(1,:); % first corner of the cell
    for ip=1:ppc(2)
      for jp=1:ppc(1)
        x(1) = x1(1) + dx*0.5 + (jp-1)*dx;
        x(2) = x1(2) + dy*0.5 + (ip-1)*dy;
        if ( inpolygon (x(1), x(2), geo(:,1), geo(:,2) ) )
          volume  = [volume;dx*dy];          
          coord   = [coord;x];
        end
      end
    end
  end
end

res.position  = coord;
res.volume    = volume;



