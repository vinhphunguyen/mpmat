function data = getCPDIPolyhedronBasis(pid,xp,particle,mesh)
% Compute CPDI-ngon (3D) shape functions and first derivatives at
% particle "pid".
%
% Inputs:
%
% pid:      particle index
% xp:       particle centroid coords
% particle: particle mesh
% mesh:     background mesh/grid
% res:      kind of global variable, 
%
% VP Nguyen, 
% 20 April 2016


% corners of the particle domain
% the name 'corners' is no longer correct as 
% we include the internal node (the particle center itself)

nodeIds = particle.element{pid};      % index
corners = particle.node(nodeIds,:);   % coordinates
faces   = particle.elem{pid}.faces;   % all faces of polyhedron particle 'pid'
cornerCount0 = size(corners,1);
numnode      = size(particle.node,1);

res.wf = zeros(numnode+1,1); res.wg = zeros(numnode+1,3);

[res,V,facesData]   = getCPDIPolyhedronData(faces,xp,particle,res);

% add centroid and face centers to the corners

corners = [corners;facesData.coord;xp];

cornerCount = size(corners,1);

wf = zeros(cornerCount,1);
wg = zeros(cornerCount,3);

for i=1:cornerCount-1 % not consider the centroid yet
    if i <= cornerCount0
        cc       = nodeIds(i);
        wf(i)    = res.wf(cc);
        wg(i,:)  = res.wg(cc,:);
    else
        cc       = i-cornerCount0;
        wf(i)    = facesData.wf(cc);
        wg(i,:)  = facesData.wg(cc,:);
    end
end

wf(end)   = res.wf(end);
wg(end,:) = res.wg(end,:);

%find elements contain the corners
elems     = zeros(cornerCount,1);    % indices of elements of particle nodes
for c=1:cornerCount
  xc        = corners(c,:);
  elems(c)  = point2ElemIndex3D(xc,mesh);
end
nodes       = unique(mesh.element(elems,:));
nodeCount   = length(nodes);
phi         = zeros(nodeCount,1);
dphi        = zeros(nodeCount,3);

% compute phi_I(xp) and first derivatives

for i=1:nodeCount
    xI = mesh.node(nodes(i),:);
    for c=1:cornerCount
        x        = corners(c,:) - xI;
        [N,~]    = getMPM3D(x,mesh.deltax,mesh.deltay,mesh.deltaz);
        phi(i)   = phi(i)    + wf(c)   * N;
        dphi(i,:)= dphi(i,:) + wg(c,:) * N;
    end    
end

det = 1/V;

data.phi  = det*phi;
data.dphi = det*dphi;
data.node = nodes;



