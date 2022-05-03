% Write PolyMesher mesh to file so that interface-gen can duplicate the
% nodes. Format is
%  NumberOfNodess 5047 NumberOfElements 720
%  NODES
%  1 0.000000  0
%  2 0.166667  0
%  ELEMENTS
%  node1 node2 node3
%  node1 node2 node3 node4
% Written by:
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia
% 25 August 2015.

function writePolyMesh(file,node,elem)

numNodes = size(node,1);
numElems = size(elem,1);

% Output files
out = fopen(file, 'wt');

%% Write headers
fprintf(out, 'NumberOfNodess %g NumberOfElements %g \n',numNodes,numElems);


%% Write node coordinates

fprintf(out, 'NODES \n');

for i=1:numNodes
    fprintf(out, '%g %f %f \n', i, node(i,1:2));
end

%% Write element data

fprintf(out, 'ELEMENTS \n');

% connectivity of every polygons
for i=1:numElems
  connect = elem{i};
  fprintf(out, '%g ', connect );
  fprintf(out, '\n');
end

% close file
fclose(out);
