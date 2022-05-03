function VTKParticlesCPDI(particles,vtuFile,data)
% Write VTK post-processing files
% for Material Point Method with CPDI2 interpolation.
% Unstructured Q4 mesh with data at element centers which are particles.
%
%
% VP Nguyen
% May 2014
% Saigon, Vietnam

node     = particles.node;
elementV = particles.elem;
etype    = particles.elemType;

dim      = size(node,2);
numNodes = size(node,1);
numCells = size(elementV,1);
x        = node;
connect  = elementV;

if iscell(connect)
  connect = cell2mat(connect);
end

% Output files

outfileVTU  = strcat(vtuFile, '.vtu');
results_vtu = fopen(outfileVTU, 'wt');

if(strcmp(etype, 'Q4') || strcmp(etype, 'Quad8') || strcmp(etype, 'Quad9'))
    numVertexesPerCell = 4;
    VTKCellCode = 9;
elseif(strcmp(etype, 'B8'))
    numVertexesPerCell = 8;
    VTKCellCode = 12;
elseif(strcmp(etype, 'T3') || strcmp(etype, 'Tri6'))
    numVertexesPerCell = 3;
    VTKCellCode = 5;
elseif(strcmp(etype, 'H4'))   % convention of Chessa's code
    numVertexesPerCell = 4;
    VTKCellCode = 10;
else
    error('Element type not known (VTKPostProcess)')
end


dof_per_vertex = 2;


%% Write headers
fprintf(results_vtu, '<VTKFile type="UnstructuredGrid"  version="0.1"   > \n');
fprintf(results_vtu, '<UnstructuredGrid> \n');
fprintf(results_vtu, '<Piece  NumberOfPoints="  %g" NumberOfCells=" %g"> \n', numNodes, numCells);

%% Write point data
fprintf(results_vtu, '<Points> \n');
if( dof_per_vertex == 1)
    fprintf(results_vtu, '<DataArray  type="Float64"  NumberOfComponents="1"  format="ascii" > \n');
else
    fprintf(results_vtu, '<DataArray  type="Float64"  NumberOfComponents="3"  format="ascii" > \n');
end

for i=1:numNodes
    if( dim == 3)
        fprintf(results_vtu, '%f ',  x(i,1:3));
    elseif(dim == 2)
        fprintf(results_vtu, '%f ',  x(i,1:2));
        fprintf(results_vtu, '0.0 ');
    end
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');
fprintf(results_vtu, '</Points> \n');

%% Print cells
fprintf(results_vtu, '<Cells> \n');

%% Print cell connectivity
fprintf(results_vtu, '<DataArray  type="Int32"  Name="connectivity"  format="ascii"> \n');

for i=1:numCells
    fprintf(results_vtu, '%g ',  connect(i,1:numVertexesPerCell)-1 );
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');

% Print cell offsets
fprintf(results_vtu, '<DataArray  type="Int32"  Name="offsets"  format="ascii"> \n');

offset = 0;
for i=1:numCells
    offset = offset + numVertexesPerCell;
    fprintf(results_vtu, '%g ', offset);
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');

% Print cell types
fprintf(results_vtu, '<DataArray  type="UInt8"  Name="types"  format="ascii"> \n');

for i=1:numCells
    fprintf(results_vtu, '%g ', VTKCellCode);
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');
fprintf(results_vtu, '</Cells> \n');

%% Print result data

sigma   = data.stress;

sigmaXX = sigma(:,1);
sigmaYY = sigma(:,2);
sigmaXY = sigma(:,3);

if size(sigma,2)==4
    sigmaVM = sigma(:,4); % von Mises stress
end


fprintf(results_vtu, '<CellData  Vectors="sigma"> \n');

if size(sigma,2)==4
    fprintf(results_vtu, '<DataArray  type="Float64"  Name="sigma" NumberOfComponents="4" format="ascii"> \n');
else
    fprintf(results_vtu, '<DataArray  type="Float64"  Name="sigma" NumberOfComponents="3" format="ascii"> \n');
end


for i=1:size(sigma,1)
    fprintf(results_vtu, '%f   ', sigmaXX(i) );
    fprintf(results_vtu, '%f   ', sigmaYY(i) );
    fprintf(results_vtu, '%f   ', sigmaXY(i) );
    
    if size(sigma,2)==4
        fprintf(results_vtu, '%f   ', sigmaVM(i) );
    end
    
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');

fprintf(results_vtu, '</CellData> \n');

% print displacement field

fprintf(results_vtu, '<CellData  Vectors="sigma"> \n');

if size(sigma,2)==4
    fprintf(results_vtu, '<DataArray  type="Float64"  Name="sigma" NumberOfComponents="4" format="ascii"> \n');
else
    fprintf(results_vtu, '<DataArray  type="Float64"  Name="sigma" NumberOfComponents="3" format="ascii"> \n');
end


for i=1:size(sigma,1)
    fprintf(results_vtu, '%f   ', sigmaXX(i) );
    fprintf(results_vtu, '%f   ', sigmaYY(i) );
    fprintf(results_vtu, '%f   ', sigmaXY(i) );
    
    if size(sigma,2)==4
        fprintf(results_vtu, '%f   ', sigmaVM(i) );
    end
    
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');

% print displacement field

if (isfield(data,'disp'))
  fprintf(results_vtu, '<DataArray  type="Float64"  Name="U" NumberOfComponents="3" format="ascii"> \n');
  for i=1:size(data.disp,1)
    fprintf(results_vtu, '%f   ', data.disp(i,:) );  
    fprintf(results_vtu, '%f   ', 0 );  
    fprintf(results_vtu, '\n');
  end
  fprintf(results_vtu, '</DataArray> \n');  
end

% print color field

if (isfield(data,'color'))
  fprintf(results_vtu, '<DataArray  type="Float64"  Name="color"  format="ascii"> \n');
  for i=1:size(data.color,1)
    fprintf(results_vtu, '%f   ', data.color(i) );  
    fprintf(results_vtu, '\n');
  end
  fprintf(results_vtu, '</DataArray> \n');  
end



fprintf(results_vtu, '</CellData> \n');

% end of VTK file

fprintf(results_vtu, '</Piece> \n');
fprintf(results_vtu, '</UnstructuredGrid> \n');
fprintf(results_vtu, '</VTKFile> \n');

fclose(results_vtu);
