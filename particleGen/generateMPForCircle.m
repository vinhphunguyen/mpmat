function [res] = generateMPForCircle(geo,ppc,mesh)
%
% Generate material points for a circle 'geo' where PPC (particle per cell)
% is given by 'ppc' and background mesh is 'mesh'.
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% August 2015.

xc   = geo.center(1);
yc   = geo.center(2);
rad  = geo.radius;

% bounding box of the circle
xmin = xc - rad;
xmax = xc + rad;
ymin = yc - rad;
ymax = yc + rad;

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

dx  = mesh.deltax/ppc(1);
dy  = mesh.deltay/ppc(2);

volume = [];
coord  = [];

ii = 1;
for i = iMin:iMax
  for j = jMin:jMax
    id  = i + mesh.numx * ( j - 1);
    res.elems(ii) = id;
    ii  = ii + 1;
    sctr = mesh.element(id,:);          %  element scatter vector
    pts  = mesh.node(sctr,:);
    x1   = pts(1,:); % first corner of the cell
    for ip=1:ppc
      for jp=1:ppc
        x(1) = x1(1) + dx*0.5 + (jp-1)*dx;
        x(2) = x1(2) + dy*0.5 + (ip-1)*dy;
        if ( (x(1) - xc)*(x(1) - xc) + (x(2) - yc)*(x(2) - yc) < rad*rad)
          volume  = [volume;dx*dy];          
          coord   = [coord;x];
        end
      end
    end
  end
end

res.position  = coord;
res.volume    = volume;


