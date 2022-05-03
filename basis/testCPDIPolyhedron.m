% test case for CPDI polyhedron
% a cube 1x1x1 with one single element

clear;

node =[0 0 0;
       1 0 0;
       1 0 1;
       0 0 1;
       0 1 0;
       1 1 0;
       1 1 1;
       0 1 1];
   

element{1} = [1 2 3 4 5 6 7 8];
elem{1}.faces{1} = [1 2 3 4];
elem{1}.faces{2} = [2 3 7 6];
elem{1}.faces{3} = [1 2 6 5];
elem{1}.faces{4} = [1 4 8 5];
elem{1}.faces{5} = [5 6 7 8];
elem{1}.faces{6} = [4 3 7 8];

numnode = size(node,1);
numelem = length(element);
pCount  = 1;

volume  = zeros(pCount,1);
coords  = zeros(pCount,3);

particles.node     = node;
particles.element  = element; % element connectivity
particles.elem     = elem;    % faces of every element

res.wf = zeros(numnode+1,1); res.wg = zeros(numnode+1,3);
for e = 1:numelem
    corners     = node(element{e},:);
    faces       = elem{e}.faces;
    coords(e,:) = centroid(corners); % center of each element=particle  
    %coords(e,:) = mean(corners); % center of each element=particle 
    [res,V]     = getCPDIPolyhedronData(faces,coords(e,:),particles,res);
   
    volume(e)   = V;  
end

ghostCell=0;
lx     = 2;
ly     = 8;
lz     = 8;
numx2  = 1;      % 2^numx2 = number of elements along X direction
numy2  = 4;      % 2^numy2 = number of elements along Y direction
numz2  = 4;      % 2^numz2 = number of elements along X direction
[mesh] = buildGrid3D(lx,ly,lz,numx2,numy2,numz2);


%res.wf = zeros(numnode,1); res.wg = zeros(numnode,3);

p=1;
data  = getCPDIPolyhedronBasis(p,coords(p,:),particles,mesh);


