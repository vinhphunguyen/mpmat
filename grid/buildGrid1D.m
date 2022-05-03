function [mesh]=buildGrid1D(L,no0, ghostCell)
% Build structured mesh for 1D domain [0,L]
% no0 = number of elements 
% ghostCell = 1: add extra cell
% ghostCell = 0: no  extra cell

no = no0;

deltax = L/no0;

if (ghostCell)
    no = no0 + 2;
    L  = L + 2*deltax;
end

% build the mesh
nodes     = linspace(0,L,no+1)';
elements  = zeros(no,2);

for ie=1:no
    elements(ie,:) = ie:ie+1;
end


% find boundaries
lNodes = 0;
rNodes = no+1;

if (ghostCell)
    lNodes=[0;1];
    rNodes=[no;no+1];
    %nodes = nodes-deltax; %normal cells within [0,L]
end

mesh.node      = nodes;
mesh.element   = elements;
mesh.deltax    = deltax;
mesh.elemCount = size(elements,1);
mesh.nodeCount = size(nodes,1);
mesh.lNodes    = lNodes;
mesh.rNodes    = rNodes;
mesh.dxInv     = 1/deltax;
mesh.ghostCell = ghostCell;

