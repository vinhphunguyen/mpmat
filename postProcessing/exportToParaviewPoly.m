function exportToParaviewPoly(xnod,glbFace,ndim,sumc,totalf,...
    ex,ey,ez,filename)

%% ***** ***** ***** Writing the vtk file ***** ***** *****
fid = fopen([filename '.vtk'],'w');

%% Writing the header from the vtk file
fprintf(fid,'# vtk DataFile Version 3.1 \n');
fprintf(fid,'This file was created by matlab source code \n');
fprintf(fid,'ASCII \n');
fprintf(fid,'DATASET POLYDATA \n');
%fprintf(fid,'MESH dimension %3.0f   Elemtype %s   Nnode %2.0f \n \n', ndim, eletyp, nno_por_elem);

%% Writing the coordinates of each node
totalNno = size(xnod,1);  % Total of nodes
fprintf(fid,'POINTS %3.0f  FLOAT \n', totalNno);
for i=1:totalNno
    fprintf(fid, [repmat('%12.5f ',1,ndim) '\n'], xnod(i,:));
end

%% Initial parameters
nef = size(glbFace,2);
cell_id = 1;  % Cell id indicate the id of the type of finite element

%% Writing the cells or nodes
totalCells = totalf;
fprintf(fid, '\n');
fprintf(fid,'POLYGONS %3.0f %3.0f \n', nef, totalCells);
for i=1:nef
    LaG = glbFace{i} ;
    nno = length(LaG);
    LaG_ = LaG - 1;
    fprintf(fid, ['%5d ' repmat('%5d ',1,nno) '\n'], nno, LaG_);
end

%% Writing the header of the stresses and deformations
nef = size(ex,1);

%% Writing the deformation Ex
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid,'POINT_DATA %3.0f \n', totalNno);
fprintf(fid,'SCALARS Ex float \n');
fprintf(fid,'LOOKUP_TABLE default \n');
for i=1:nef
    fprintf(fid,'%12.5g \n', ex(i));
end

fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid,'POINT_DATA %3.0f \n', totalNno);
if ~isempty(ey)
%% Writing the deformation Ey
fprintf(fid,'SCALARS Ey float \n');
fprintf(fid,'LOOKUP_TABLE default \n');
for i=1:nef
    fprintf(fid,'%12.5g \n', ey(i));
end
end

fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid,'POINT_DATA %3.0f \n', totalNno);
if ~isempty(ez)
%% Writing the deformation Ez
fprintf(fid, '\n');
fprintf(fid,'SCALARS Ez float \n');
fprintf(fid,'LOOKUP_TABLE default \n');
for i=1:nef
    fprintf(fid,'%12.5g \n', ez(i));
end
end

fclose(fid);
return;
