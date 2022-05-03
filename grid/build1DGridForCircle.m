function [mesh]=build1DGridForCircle(center, radius,elemCount)
% build a mesh of two-noded bar elements for a circle
% defined by center and radius.
% Phu Nguyen
% Monash University
% 21 July 2016

% first circle at (0,0)
alpha    = 2*pi/elemCount; 
nodes    = zeros(elemCount,2);
elements = zeros(elemCount,2);

for i=1:elemCount
    nodes(i,:)    = [radius*cos((i-1)*alpha) radius*sin((i-1)*alpha)];
    elements(i,:) = [i i+1]; 
end

elements(elemCount,2) = 1;

% move nodes to correct position
nodes(:,1) = nodes(:,1) + center(1);
nodes(:,2) = nodes(:,2) + center(2);

% write output
mesh.node      = nodes;
mesh.element   = elements;
mesh.elemCount = size(elements,1);
mesh.nodeCount = size(nodes,1);
mesh.type      = 'L2'; % to use Chessa;s plot

% to plot the mesh 
% plot_mesh(mesh.node,mesh.element,mesh.type,'r-*',2)
