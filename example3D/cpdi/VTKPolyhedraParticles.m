% Write Material Point Method results to VTK files
% which can be processed by ViSIT and/or Paraview.
% This vtp file is for polyhedra CPDI MPM simulations.
% Written by:
% Vinh Phu Nguyen
% Monash University
% 21 April 2016.

function VTKPolyhedraParticles(particles,vtpFile,data)

% particles:  particle information
% vtpFile:    VTK file to be written
% data:       data to be written to file

node     = particles.node;
elem     = particles.elem;

numNodes = size(node,1);
numElems = length(elem);
numFaces = length(particles.face.vertex_indices);

% Output files

outfileVTU  = strcat(vtpFile, '.ply');
results_vtu = fopen(outfileVTU, 'wt');

%% Write headers
fprintf(results_vtu, 'ply \n');
fprintf(results_vtu, 'format ascii 1.0 \n');
fprintf(results_vtu, 'element vertex  %g\n', numNodes);
fprintf(results_vtu, 'property float x \n');
fprintf(results_vtu, 'property float y \n');
fprintf(results_vtu, 'property float z \n');
fprintf(results_vtu, 'element face   %g\n', numFaces);
fprintf(results_vtu, 'property list uchar int vertex_indices \n');
fprintf(results_vtu, 'element cell   %g\n', numElems);
fprintf(results_vtu, 'property list uchar int face_indices \n');
fprintf(results_vtu, 'end_header \n');
%% Write point(vertices) coordinates

for i=1:numNodes
    fprintf(results_vtu, '%f %f %f \n',  node(i,1:3));
end

%% Write faces

for f=1:numFaces
    face = particles.face.vertex_indices{f};
    fprintf(results_vtu, '%g %g ', length(face), face );
    fprintf(results_vtu, '\n');
end



% close file
fclose(results_vtu);
