function data = cpdi2Triangle(pid,particle,mesh)
% Compute CPDI2 (Triangle) shape functions and first derivatives at
% particle "pid".
%
% Inputs:
%
% pid:      particle index
% particle: particle mesh
% mesh:     background mesh/grid
%
% VP Nguyen
% May, 2014
% Saigon, Vietnam

% three corners of the particle domain

nodeIds = particle.elem(pid,:);
corners = particle.node(nodeIds,:);

% function and gradient weights

y13   = (corners(1,2)-corners(3,2));
y23   = (corners(2,2)-corners(3,2));
y12   = (corners(1,2)-corners(2,2));
x12   = (corners(1,1)-corners(2,1));
x31   = (corners(3,1)-corners(1,1));
x32   = (corners(3,1)-corners(2,1));

wg    = zeros(3,2);
wf    = [1/3 1/3 1/3];

wg(1,:) = 0.5*[ y23  x32];
wg(2,:) = 0.5*[-y13 -x31];
wg(3,:) = 0.5*[ y12 -x12];

% particle domain area

Vp     = 0.5*(y13*x12+y12*x31);

elems  = zeros(3,1); % indices of elements of 3 corners

% find elements contain the corners

for c=1:3
    xc        = corners(c,:);
    elems(c)  = point2ElemIndex(xc,mesh);
end

% nodes I where phi_I(xp) are non-zero

nodes = unique(mesh.element(elems,:));

% compute phi_I(xp) and first derivatives

nodeCount = length(nodes);
phi       = zeros(nodeCount,1);
dphi      = zeros(nodeCount,2);

for i=1:nodeCount
    xI = mesh.node(nodes(i),:);
    for c=1:3
        x        = corners(c,:) - xI;
        [N,dNdx] = getMPM2D(x,mesh.deltax,mesh.deltay);
        phi(i)   = phi(i)    + wf(c)  *N;
        dphi(i,:)= dphi(i,:) + wg(c,:)*N;
    end    
end

dphi      = dphi/(Vp);

data.phi  = phi;
data.dphi = dphi;
data.node = nodes;



