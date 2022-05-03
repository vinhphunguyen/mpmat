function [res] = generateMPQuadDiffCircles(quad,circle1,circle2,ppc,mesh)
%
% Generate material points for a quadrilateral 'quad' not within a circle 
% where PPC (particle per cell)
% is given by 'ppc' and background mesh is 'mesh'.
% quad is a matrix 4x2.
%
% Use built-in function 'inpolygon' to check if a point is within a
% polygon.
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% August 2015.

% bounding box of the quad
mi = min(quad);
ma = max(quad);

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

dx  = mesh.deltax/(ppc);
dy  = dx;

volume = [];
coord  = [];

xc1  = circle1.center(1);
yc1  = circle1.center(2);
rc12 = circle1.radius * circle1.radius;


xc2  = circle2.center(1);
yc2  = circle2.center(2);
rc22 = circle2.radius * circle2.radius;

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
        if ( inpolygon (x(1), x(2), quad(:,1), quad(:,2) ) == 0 )          
          continue
        elseif ( (x(1)-xc1)*(x(1)-xc1) + (x(2)-yc1)*(x(2)-yc1) > rc12 ) && ...
               ( (x(1)-xc2)*(x(1)-xc2) + (x(2)-yc2)*(x(2)-yc2) > rc22 ) 
          volume  = [volume;dx*dy];          
          coord   = [coord;x];
        end
      end
    end
  end
end

res.position  = coord;
res.volume    = volume;



