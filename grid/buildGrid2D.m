function [mesh]=buildGrid2D(Lx,Ly,noX0,noY0, ghostCell)
% Build structured mesh for rectangle (Lx,Ly)
% noX0 = number of elements along x direction
% noY0 = number of elements along y direction
% ghostCell = 1: add extra cell
% ghostCell = 0: no  extra cell
%
% Vinh Phu Nguyen
% June 2014, Saigon, Vietnam.

noX = noX0;
noY = noY0;

deltax = Lx/noX0;
deltay = Ly/noY0;

if (ghostCell)
    noX = noX0 + 2;
    noY = noY0 + 2;
    Lx  = Lx + 2*deltax;
    Ly  = Ly + 2*deltay;
end

% build the mesh
nnx  = noX+1;
nny  = noY+1;
node = square_node_array([0 0],[Lx 0],[Lx Ly],[0 Ly],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element =make_elem(node_pattern,noX,noY,inc_u,inc_v);

% find boundaries
eps = 1e-12;
lNodes = find(abs(node(:,1))    < eps);
rNodes = find(abs(node(:,1)-Lx) < eps);
tNodes = find(abs(node(:,2)-Ly) < eps);
bNodes = find(abs(node(:,2))    < eps);

if (ghostCell)
    lNodes0 = find(abs(node(:,1)-deltax)    < eps);
    rNodes0 = find(abs(node(:,1)-Lx+deltax) < eps);
    tNodes0 = find(abs(node(:,2)-Ly+deltay) < eps);
    bNodes0 = find(abs(node(:,2)-deltay)    < eps);
    
    lNodes=[lNodes;lNodes0];
    rNodes=[rNodes;rNodes0];
    tNodes=[tNodes;tNodes0];
    bNodes=[bNodes;bNodes0];
end

% build nodal support

nodalSup  = cell(size(node,1),1);

for e=1:size(element,1)
  sctr = element(e,:);
  for i=1:4
    nodei = sctr(i);
    nodalSup{nodei}=[nodalSup{nodei}; e];
  end
end

mesh.node      = node;
mesh.element   = element;
mesh.deltax    = deltax;
mesh.deltay    = deltay;
mesh.elemCount = size(element,1);
mesh.nodeCount = size(node,1);
mesh.numx      = noX;
mesh.numy      = noY;
mesh.lNodes    = lNodes;
mesh.rNodes    = rNodes;
mesh.tNodes    = tNodes;
mesh.bNodes    = bNodes;
mesh.dxInv     = 1/deltax;
mesh.dyInv     = 1/deltay;
mesh.nodalSup  = nodalSup;

