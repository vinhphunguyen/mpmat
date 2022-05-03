function data = getCPDIPolygonBasis(pid,input,particle,mesh)
% Compute CPDI-ngon shape functions and first derivatives at
% particle "pid".
%
% Inputs:
%
% pid:      particle index
% particle: particle mesh
% mesh:     background mesh/grid
%
% VP Nguyen, The University of Adelaide, Australia.
% August, 2015


% corners of the particle domain
% the name 'corners' is no longer correct as 
% we include the internal node (the particle center itself)

nodeIds = particle.elem{pid};         % index
corners = particle.node(nodeIds,:);   % coordinates
xp      = mean(corners);              % particle xp

% add the centroid of polygon to the corners!!!
corners = [corners;xp];

%Vp     = input.Vp;
nodes  = input.nodes;
wf     = input.wf;
wg     = input.wg;

% compute phi_I(xp) and first derivatives

cornerCount = size(corners,1);
nodeCount   = length(nodes);
phi         = zeros(nodeCount,1);
dphi        = zeros(nodeCount,2);

for i=1:nodeCount
    xI = mesh.node(nodes(i),:);
    for c=1:cornerCount
        x        = corners(c,:) - xI;
        [N,~]    = getMPM2D(x,mesh.deltax,mesh.deltay);
        phi(i)   = phi(i)    + wf(c)  *N;
        dphi(i,:)= dphi(i,:) + wg(c,:)*N;
    end    
end

data.phi  = phi;
data.dphi = dphi;
data.node = nodes;



