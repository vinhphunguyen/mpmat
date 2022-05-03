function data = getCPDIQuadBasis(pid,input,particle,mesh)
% Compute CPDI2 shape functions and first derivatives at
% particle "pid".
%
% Inputs:
%
% pid:      particle index
% particle: particle mesh
% mesh:     background mesh/grid
%
% VP Nguyen
% July, 2014
% Saigon, Vietnam

% four corners of the particle domain

nodeIds = particle.elem(pid,:);
corners = particle.node(nodeIds,:);

Vp     = input.Vp;
nodes  = input.nodes;
wf     = input.wf;
wg     = input.wg;

% compute phi_I(xp) and first derivatives

nodeCount = length(nodes);
phi       = zeros(nodeCount,1);
dphi      = zeros(nodeCount,2);

for i=1:nodeCount
    xI = mesh.node(nodes(i),:);
    for c=1:4
        x        = corners(c,:) - xI;
        [N,dNdx] = getMPM2D(x,mesh.deltax,mesh.deltay);
        phi(i)   = phi(i)    + wf(c)  *N;
        dphi(i,:)= dphi(i,:) + wg(c,:)*N;
    end    
end

% done at particle data
%phi       = phi /(36*Vp);
%dphi      = dphi/(2*Vp);

data.phi  = phi;
data.dphi = dphi;
data.node = nodes;



