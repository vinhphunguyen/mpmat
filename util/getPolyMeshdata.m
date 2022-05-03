function [node,element,Elements] = getPolyMeshdata(fileName)

% purpose: to get the polyhedron mesh data. the data is stored in a ply
% format. Open source software NEPER (http://neper.sourceforge.net) is used
% to generate the polyhedron mesh
% Ref: Large scale 3D random polycrystals for the finite element method:
% generation, meshing and remeshing, CMAME,v200, 1729--1745, 2011
%
% Input:
%   fileName - name of the file which has the mesh information. The
%   extension is ply
%
% Outputs:
%   node - nodal coordinates
%   element - element connectivity
%   Elements - a data structure containing the following
%       Elements.vertex.i (i=x,y,z) - nodal coordinates
%       Elements.face.vertex_indices - vertex nodes that build up the 
%       face, i.e.,face connectivity
%       Elements.cell.face_indices - face index that build up the element
%
% Sundar, UNSW, 2014
%--------------------------------------------------------------------------

[Elements,varargout] = plyread(fileName);

tote = size(Elements.cell.face_indices,1);

% now get the element connectivity...this is done by looping over the
% number of elements. Elements.cell.face_indices has the index of the faces
% from which the element connectivity needs to be done...
for iel = 1:tote
    
    % get the current face indices of the current element...remember the
    % index obtained from ply file starts with 0, add 1 to the vector so
    % that matlab does not throw an error
    CurFaceIndex = Elements.cell.face_indices{iel} + 1;
    
    % the variable CurFaceIndex has the indices of the faces. From this,
    % get the vertex numbers. This is stored in
    % Elements.face.vertex_indices. Remember the index obtained from the
    % ply file starts with 0, add 1 to the vector so that matlab does not
    % throw an error. Also need to keep track is whether a particular
    % vertex number is already in the element connectivity, if so do not
    % add again. Loop over the Face indices to construct the element
    % connectivity
    
    tempE = [];
    for in = 1:length(CurFaceIndex)
        
        % get the first index.
        CFindex = CurFaceIndex(in);
        
        % get the vertex number of the current face
        vertexN = Elements.face.vertex_indices{CFindex} + 1;
        
        if in == 1
            tempE = [tempE vertexN];
        end
        
        if in > 1
            
            % now loop over this and add to the element connectivity. keep
            % track if a particular verex number is already assigned.
            for im = 1:length(vertexN)
                
                % check if the current vertex is already in the list
                % maintained in tempE. if it is not then add to the list.
                % else do nothing
                if ~ismember(vertexN(im),tempE)
                    tempE = [tempE vertexN(im)];
                end
            end
        end
    end
    
    element{iel} = tempE;
end


node = [Elements.vertex.x Elements.vertex.y Elements.vertex.z];

end