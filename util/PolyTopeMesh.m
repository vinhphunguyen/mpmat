clear all
close all
clc

% load the mesh. the mesh is generated with neper open source software. the
% output of which is written in ply format. read the ply format to get the
% nodal coordinates, element and face connectivity.
fprintf(1,'....Load the mesh\n');
fileName = 'temp1.ply';

[node,element,Elements] = getPolyMeshdata(fileName);
numelem = size(element,2)
numnode = size(node,1)

%... vertex indices for all the polytopes...
Faces = Elements.face.vertex_indices;

% loop over elements...
for iel = 1:numelem
    
    % face ids of the current element....
    faceid = Elements.cell.face_indices{iel} + 1;
    
    clear faces
    nds = [];
    % get the faces connecting this particular element
    for in = 1:length(faceid)
        
        % corner indexes of the current polytope
        tm = Faces{faceid(in)} + 1;
        
        % nodal coordinates of the current face
        nds = [nds; node(tm,:)] ;
        
        faces{in} = tm;
    end
    
    %...this saves the list of corner indexes for the current polytope
    elem{iel,1}.faces = faces;
end
