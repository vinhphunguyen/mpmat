function data = cpdi2D(xp,r1p,r2p,mesh)
% Compute CPDI shape functions and first derivatives at
% particle xp.
%
% Inputs:
%
% r1p, r2p: particle domain vectors
% xp:       particle position
% mesh:     background mesh/grid
%
% VP Nguyen
% May, 2014
% Saigon, Vietnam

% four corners of the particle domain

x1 = xp - r1p - r2p;
x2 = xp + r1p - r2p;
x3 = xp + r1p + r2p;
x4 = xp - r1p + r2p;

% gradient weights

w    = zeros(4,2);

w(1,:) = [r1p(2)-r2p(2) r2p(1)-r1p(1)];
w(2,:) = [r1p(2)+r2p(2) -r1p(1)-r2p(1)];
w(3,:) = -w(1,:);
w(4,:) = -w(2,:);

% particle domain area

Vp     = 4*norm(cross([r1p 0],[r2p 0]));

data.x = [x1;x2;x3;x4];
elems  = zeros(4,1); % indices of elements of 4 corners

% find elements contain the corners

for c=1:4
    xc        = data.x(c,:);
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
    for c=1:4
        x        = data.x(c,:) - xI;
        [N,dNdx] = getMPM2D(x,mesh.deltax,mesh.deltay);
        phi(i)   = phi(i)    + N;
        dphi(i,:)= dphi(i,:) + w(c,:)*N;
    end    
end

phi       = phi *0.25;
dphi      = dphi/Vp;

data.phi  = phi;
data.dphi = dphi;
data.node = nodes;



