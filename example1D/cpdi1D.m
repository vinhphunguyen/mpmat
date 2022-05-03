function data = cpdi1D(xp,rp,mesh)
% Compute CPDI shape functions and first derivatives at
% particle xp.
%
% Inputs:
%
% rp:       particle domain vector
% xp:       particle position
% mesh:     background mesh/grid
%
% VP Nguyen
% May, 2014
% Saigon, Vietnam

% two ends of the particle domain

x1 = xp - rp;
x2 = xp + rp;

% gradient weights

w    = zeros(2,1);

w(1) = r1p(2)-r2p(2);
w(2) = r1p(2)+r2p(2);


% particle domain area

Vp     = rp;

data.x = [x1;x2];
elems  = zeros(2,1); % indices of elements of 4 corners

% find elements contain the corners

for c=1:2
    xc        = data.x(c,:);
    elems(c)  = point2ElemIndex(xc,mesh);
end

% nodes I where phi_I(xp) are non-zero

nodes = unique(mesh.element(elems,:));

% compute phi_I(xp) and first derivatives

nodeCount = length(nodes);
phi       = zeros(nodeCount,1);
dphi      = zeros(nodeCount,1);

for i=1:nodeCount
    xI = mesh.node(nodes(i),:);
    for c=1:2
        x        = data.x(c,:) - xI;
        [N,dNdx] = getMPM(x,mesh.deltax);
        phi(i)   = phi(i)    + N;
        dphi(i)  = dphi(i)   + w(c)*N;
    end    
end

phi       = phi *0.5;
dphi      = dphi*(0.5/Vp);

data.phi  = phi;
data.dphi = dphi;
data.node = nodes;



