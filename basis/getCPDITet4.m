function data = getCPDITet4(pid,particle,mesh)
% Compute CPDI2 (linear tetrahedron) shape functions and first derivatives at
% particle "pid".
%
% Inputs:
%
% pid:      particle index
% particle: particle mesh
% mesh:     background mesh/grid
%
% VP Nguyen
% 13 April, 2016
% Monash University
% Clayton, Australia

% four corners of the particle domain

nodeIds = particle.elem(pid,:);
corners = particle.node(nodeIds,:);

% function and gradient weights

x1   = corners(1,1); y1   = corners(1,2); z1   = corners(1,3);
x2   = corners(2,1); y2   = corners(2,2); z2   = corners(2,3);
x3   = corners(3,1); y3   = corners(3,2); z3   = corners(3,3);
x4   = corners(4,1); y4   = corners(4,2); z4   = corners(4,3);

x21 = x2 - x1; x32 = x3 - x2; x43 = x4 - x3; x42 = x4 - x2; x13 = x1 - x3; x31 = -x13; x34=-x43; x14=x1-x4; x24=-x42; x12=-x21;
y23 = y2 - y3; y34 = y3 - y4; y12 = y1 - y2; y42 = y4 - y2; y32 = -y23; y31 = y3 - y1; y43=-y34; y13=-y31; y24=-y42; y14=y1-y4; y21=-y12;
z34 = z3 - z4; z23 = z2 - z3; z12 = z1 - z2; z42 = z4 - z2; z32 = -z23; z43 = -z34; z13 = z1 - z3; z31 = -z13; z14=z1-z4;z24=-z42; z21=-z12;

V   = x21*(y23*z34-y34*z23) + x32*(y34*z12-y12*z34) + x43*(y12*z23-y23*z12);
a1  = y42*z32 - y32*z42; b1 = x32*z42 - x42*z32; c1 = x42*y32 - x32*y42;
a2  = y31*z43 - y34*z13; b2 = x43*z31 - x13*z34; c2 = x31*y43 - x34*y13;
a3  = y24*z14 - y14*z24; b3 = x14*z24 - x24*z14; c3 = x24*y14 - x14*y24;
a4  = y13*z21 - y12*z31; b4 = x21*z13 - x31*z12; c4 = x13*y21 - x12*y31;

wg    = zeros(4,3);
wf    = [1/4 1/4 1/4 1/4];

det   = 1/V; % 1/(6V)

wg(1,:) = det*[a1 b1 c1];
wg(2,:) = det*[a2 b2 c2];
wg(3,:) = det*[a3 b3 c3];
wg(4,:) = det*[a4 b4 c4];

% particle domain area

elems  = zeros(4,1); % indices of elements of 4 corners

% find elements contain the corners

for c=1:4
    xc        = corners(c,:);
    elems(c)  = point2ElemIndex3D(xc,mesh);
end

% nodes I where phi_I(xp) are non-zero

nodes = unique(mesh.element(elems,:));

% compute phi_I(xp) and first derivatives

nodeCount = length(nodes);
phi       = zeros(nodeCount,1);
dphi      = zeros(nodeCount,3);

for i=1:nodeCount
    xI = mesh.node(nodes(i),:);
    for c=1:4
        x        = corners(c,:) - xI;
        [N,~]    = getMPM3D(x,mesh.deltax,mesh.deltay,mesh.deltaz);
        phi(i)   = phi(i)    + wf(c)  *N;
        dphi(i,:)= dphi(i,:) + wg(c,:)*N;
    end    
end

%dphi      = dphi/(Vp);

data.phi  = phi;
data.dphi = dphi;
data.node = nodes;



