function data = getCPDILineData(pid,particle,mesh)
% Compute CPDI-L2 (particle domain is a two-node line element) data at particle "pid".
%
% Inputs:
%
% pid:      particle index
% particle: particle mesh
% mesh:     background mesh/grid
%
% Output:
%
% data.wf:    function weights
% data.wf:    gradient weights
% data.nodes: indices of nodes influencing particle "pid"
%
% VP Nguyen
% September, 2015
% The University of Adelaide, Australia.

% two corners of the particle domain

nodeIds = particle.elem(pid,:);
corners = particle.node(nodeIds,:);

% particle domain area

Vp     = norm(corners(1,:) - corners(2,:));

% function and gradient weights
wg    = zeros(2,2);
wf    = [1/2 1/2];

wg(1,:) = [-1 1];
wg(2,:) = [corners(3,2)-corners(1,2) corners(1,1)-corners(3,1)];


wg(:,:) = (1/(2*Vp))*wg(:,:);

elems  = zeros(4,1); % indices of elements of 4 corners

% find elements contain the corners

for c=1:4
    xc        = corners(c,:);
    elems(c)  = point2ElemIndex(xc,mesh);
end

% nodes I where phi_I(xp) are non-zero

nodes = unique(mesh.element(elems,:));

data.nodes = nodes;
data.wf    = wf;
data.wg    = wg;
data.Vp    = Vp;



