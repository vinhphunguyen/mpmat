function data = getCPDIPolygonData(pid,particle,mesh)
% Compute CPDI-ngon data at particle "pid".
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
% 22 August, 2015
% The University of Adelaide

% corners of the particle domain

nodeIds   = particle.elem{pid};
corners   = particle.node(nodeIds,:);% coordinates of corners (vertices)
nodeCount = length(nodeIds);
xp        = mean(corners);           % centroid of the particle domain
% Note that the centroid is considered as a corner of the sub-triangles
% Therefore, there are (nodeCount+1) weights, and the centroid part is put
% at the end.
wf        = zeros(nodeCount+1,1);    % function weights
wg        = zeros(nodeCount+1,2);    % gradient weights
elems     = zeros(nodeCount+1,1);    % indices of elements of particle nodes

A = 0;

for i = 1:nodeCount
  j    = rem(i,nodeCount) + 1;
  x1   = corners(i,1); y1   = corners(i,2);
  x2   = corners(j,1); y2   = corners(j,2);
  x3   = xp(1);        y3   = xp(2);
  x21  = x2 - x1; 
  y31  = y3 - y1;
  y21  = y2 - y1;
  x31  = x3 - x1;
  area = 0.5*( x21*y31 - y21*x31 );  % area of sub-triangle
  
  A    = A + area; 
  
  area3 = area/3;
  
  wf(i)   = wf(i)   + area3;
  wf(j)   = wf(j)   + area3;
  wf(end) = wf(end) + area3;
  
  wg(i,1) = wg(i,1) + (0.5)*(y2-y3); % x-derivative
  wg(i,2) = wg(i,2) + (0.5)*(x3-x2); % y-derivative
  
  wg(j,1) = wg(j,1) + (0.5)*(y3-y1);
  wg(j,2) = wg(j,2) + (0.5)*(x1-x3);
  
  wg(end,1) = wg(end,1) + (0.5)*(y1-y2);
  wg(end,2) = wg(end,2) + (0.5)*(x2-x1);
  
  elems(i)  = point2ElemIndex([x1 y1],mesh);
end

wf = wf/A;
wg = wg/A;

elems(end)  = point2ElemIndex(xp,mesh);
nodes       = unique(mesh.element(elems,:));

data.nodes = nodes;
data.wf    = wf;
data.wg    = wg;
data.Vp    = A;

