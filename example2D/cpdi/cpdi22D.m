function data = cpdi22D(pid,particle,mesh)
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
% May, 2014
% Saigon, Vietnam

% four corners of the particle domain

nodeIds = particle.elem(pid,:);
corners = particle.node(nodeIds,:);

% particle domain area

Vp     = 0.5*   ( corners(1,1)*corners(2,2)  - corners(2,1)*corners(1,2) ...
                + corners(2,1)*corners(3,2)  - corners(3,1)*corners(2,2) ...
                + corners(3,1)*corners(4,2)  - corners(4,1)*corners(3,2) ...
                + corners(4,1)*corners(1,2)  - corners(1,1)*corners(4,2) );

% function and gradient weights

c1   = (corners(2,1)-corners(1,1))*(corners(4,2)-corners(1,2)) - ...
       (corners(2,2)-corners(1,2))*(corners(4,1)-corners(1,1));
   
c2   = (corners(2,1)-corners(1,1))*(corners(3,2)-corners(2,2)) - ...
       (corners(2,2)-corners(1,2))*(corners(3,1)-corners(2,1));
   
c3   = (corners(3,1)-corners(4,1))*(corners(4,2)-corners(1,2)) - ...
       (corners(3,2)-corners(4,2))*(corners(4,1)-corners(1,1));
   
c4   = (corners(3,1)-corners(4,1))*(corners(3,2)-corners(2,2)) - ...
       (corners(3,2)-corners(4,2))*(corners(3,1)-corners(2,1));   

% a    = (corners(4,1)-corners(1,1))*(corners(2,2)-corners(3,2)) - ...
%        (corners(2,1)-corners(3,1))*(corners(4,2)-corners(1,2));
% 
% b    = (corners(3,1)-corners(4,1))*(corners(1,2)-corners(2,2)) - ...
%        (corners(1,1)-corners(2,1))*(corners(3,2)-corners(4,2));


wg    = zeros(4,2);

%wf    = [6*Vp-a-b 6*Vp-a+b 6*Vp+a+b 6*Vp+a-b];

wf    = [4*c1+2*c2+2*c3+c4 2*c1+4*c2+c3+2*c4 c1+2*c2+2*c3+4*c4 2*c1+c2+4*c3+2*c4];

wg(1,:) = [corners(2,2)-corners(4,2) corners(4,1)-corners(2,1)];
wg(2,:) = [corners(3,2)-corners(1,2) corners(1,1)-corners(3,1)];
wg(3,:) = -wg(1,:);
wg(4,:) = -wg(2,:);



elems  = zeros(4,1); % indices of elements of 4 corners

% find elements contain the corners

for c=1:4
    xc        = corners(c,:);
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
        x        = corners(c,:) - xI;
        [N,dNdx] = getMPM2D(x,mesh.deltax,mesh.deltay);
        phi(i)   = phi(i)    + wf(c)  *N;
        dphi(i,:)= dphi(i,:) + wg(c,:)*N;
    end    
end

%phi       = phi /(24*Vp);
phi       = phi /(36*Vp);
dphi      = dphi/(2*Vp);

data.phi  = phi;
data.dphi = dphi;
data.node = nodes;



