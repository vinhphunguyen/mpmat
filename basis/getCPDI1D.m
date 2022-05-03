function data = getCPDI1D(pid,particle,mesh)
% Compute CPDI2 shape functions and first derivatives at
% particle "pid" for two-node particle domains.
%
% Inputs:
%
% pid:      particle index
% particle: particle mesh
% mesh:     background mesh/grid
%
% VP Nguyen
% July, 2016
% Monash University

% two corners of the particle domain

nodeIds = particle.element(pid,:);
corners = particle.node(nodeIds);

% particle domain area

Vp     = corners(2)-corners(1);

% function and gradient weights

wf    = [0.5 0.5];
wg    = [-1/Vp 1/Vp];


elems  = zeros(2,1); % indices of elements of two corners

% find elements contain the corners

for c=1:2
    xc        = corners(c);
    elems(c)  = point2ElemIndex1D(xc,mesh);
end

% nodes I where phi_I(xp) are non-zero
elems;
nodes = unique(mesh.element(elems,:));

% compute phi_I(xp) and first derivatives

nodeCount = length(nodes);
phi       = zeros(nodeCount,1);
dphi      = zeros(nodeCount,1);

for i=1:nodeCount
    xI = mesh.node(nodes(i));
    for c=1:2
        x        = corners(c) - xI;
        [N,~]    = getMPM(x,mesh.deltax);
        phi(i)   = phi(i)  + wf(c)*N;
        dphi(i)  = dphi(i) + wg(c)*N;
    end    
end

data.phi  = phi;
data.dphi = dphi;
data.node = nodes;



