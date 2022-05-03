% Write Material Point Method results to VTK files
% which can be processed by ViSIT and/or Paraview.
% This vtp file is for polygonal CPDI MPM simulations.
% Written by:
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia
% 22 August 2015.

function VTKPolygonParticles(particles,vtpFile,data)

% particles:  particle information
% vtpFile:    VTK file to be written
% data:       data to be written to file

node     = particles.node;
elem     = particles.elem;

dim      = size(node,2);
numNodes = size(node,1);
numElems = size(elem,1);

x        = node;
if (dim==2)
    x(:,3) = 0;
end

% Output files

outfileVTU  = strcat(vtpFile, '.vtp');
results_vtu = fopen(outfileVTU, 'wt');

%% Write headers
fprintf(results_vtu, '<VTKFile type="PolyData"  version="0.1"   > \n');
fprintf(results_vtu, '<PolyData> \n');
fprintf(results_vtu, '<Piece  NumberOfPoints="  %g" NumberOfVerts=" %g" NumberOfLines=" %g" NumberOfStrips=" %g" NumberOfPolys=" %g"> \n',...
    numNodes, 0,0,0,numElems);


%% Write point coordinates

fprintf(results_vtu, '<Points> \n');
fprintf(results_vtu, '<DataArray  type="Float64"  NumberOfComponents="3"  format="ascii" > \n');

for i=1:numNodes
    fprintf(results_vtu, '%f %f %f \n',  x(i,1:3));
end

fprintf(results_vtu, '</DataArray> \n');
fprintf(results_vtu, '</Points> \n');

%% Write poly data

fprintf(results_vtu, '<Polys> \n');

% connectivity of every polygons
fprintf(results_vtu, '<DataArray  type="Int32"  Name="connectivity"  format="ascii"> \n');
for i=1:numElems
  connect = elem{i}-1;
  fprintf(results_vtu, '%g ', connect );
  fprintf(results_vtu, '\n');
end
fprintf(results_vtu, '</DataArray> \n');

% Print cell offsets
fprintf(results_vtu, '<DataArray  type="Int32"  Name="offsets"  format="ascii"> \n');

offset = 0;
for i=1:numElems      
    numVertexesPerCell = length(elem{i});
    offset = offset + numVertexesPerCell;
    fprintf(results_vtu, '%g ', offset);
    fprintf(results_vtu, '\n');
end
fprintf(results_vtu, '</DataArray> \n');

fprintf(results_vtu, '</Polys> \n');

%% write point data

% fprintf(results_vtu, '<CellData  Vectors="U"> \n');
% 
% if (isfield(data,'disp'))
%   fprintf(results_vtu, '<DataArray  type="Float64"  Name="U" NumberOfComponents="3" format="ascii"> \n');
%   for i=1:size(data.disp,1)
%     fprintf(results_vtu, '%f   ', data.disp(i,:) );  
%     fprintf(results_vtu, '%f   ', 0              );  
%     fprintf(results_vtu, '\n');
%   end
%   fprintf(results_vtu, '</DataArray> \n');
% end
% 
% fprintf(results_vtu, '</CellData> \n');
% 
% fprintf(results_vtu, '<CellData  Vectors="color"> \n');
% 
% if (isfield(data,'color'))
%   fprintf(results_vtu, '<DataArray  type="Float64"  Name="color"  format="ascii"> \n');
%   for i=1:size(data.color,1)
%     fprintf(results_vtu, '%f   ', data.color(i) );  
%     fprintf(results_vtu, '\n');
%   end
%   fprintf(results_vtu, '</DataArray> \n');
% end
% 
% fprintf(results_vtu, '</CellData> \n');

%% end of VTK file

fprintf(results_vtu, '</Piece> \n');
fprintf(results_vtu, '</PolyData> \n');
fprintf(results_vtu, '</VTKFile> \n');

% close file
fclose(results_vtu);
