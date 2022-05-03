%
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% Particle domains are represented by polygons obtained from PolyMesher.
%
% Cantilever beam with large deformation.
%
% Vinh Phu Nguyen
% Monash University
% 13 April 2016.

%%

addpath ../../nurbs/nurbs-geopdes/inst
addpath ../../nurbs/nurbs-util/
addpath ../../util/
addpath ../../fem/
addpath ../../postProcessing/
addpath ../../constitutiveModels/
addpath ../../grid/
addpath ../../basis/
addpath ../../PolyMesher/
addpath ../../geoMesh/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
%
E   = 1e6;           % Young's modulus
nu  = 0.3;           % Poisson ratio
rho = 1050;          % density
K   = E/3/(1-2*nu);    % bulk modulus
mu    = E/2/(1+nu);% shear modulus
lambda = K - 2/3*mu;

g     = 10; % gravity
bodyf = [0 -g];

I  = [1 0;0 1];

interval     = 1;
vtkFileName  = 'cantileverPoly';
vtkFileName1 = '../../results/cpdi/cantileverPolyGrid';

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell=0;
lx     = 8;
ly     = 8;
numx2  = 16;      % number of elements along X direction
numy2  = 16;      % number of elements along Y direction
[mesh] = buildGrid2D(lx,ly,numx2,numy2, ghostCell);
element= mesh.element;
node   = mesh.node;

% find boundary nodes

fixNodes=find(abs(node(:,1)-mesh.deltax)<1e-10);
%%   particle distribution
% which is in this case a rectangle
NElem = 120;
[Node,Element,Supp,Load,P]=PolyMesher(@CantileverDomain,NElem,80);

% find marked corner (recording displacement)
markedNode1 = find(abs(Node(:,1)-4)<1e-6);
markedNode2 = find(abs(Node(:,2)-0)<1e-6);

markedNode  = intersect(markedNode1,markedNode2);

Node(:,1) = Node(:,1) +  mesh.deltax;
Node(:,2) = Node(:,2) + 3.5;


markedYCoord = Node(markedNode,2);

% store the particle mesh into a structure for convenience
particles.node     = Node;
particles.elem     = Element;

pCount  = size(Element,1);                % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
coords  = zeros(pCount,2);
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stressx
strain  = zeros(pCount,3);                % strain
velo    = zeros(pCount,2);
position= zeros(pCount,2);                % particle position
color    = zeros(pCount,1);
% CPDI data
nodeid  = cell(pCount,1);                 % nodes affect particle 'p'
funcW   = cell(pCount,1);                 % function weights of 'p'
gradW   = cell(pCount,2);                 % gradient weights of 'p', 
volum   = cell(pCount,1);                 % store particle domain volume

% particle mass, volume and initial velocity
for e = 1:NElem
  vx  = Node(Element{e},1); 
  vy  = Node(Element{e},2); 
  nv  = length(Element{e});
  vxS = vx([2:nv 1]); 
  vyS = vy([2:nv 1]); %Shifted vertices
  temp        = vx.*vyS - vy.*vxS;
  a           = 0.5*sum(temp);
  volume(e)   = a;
  mass(e)     = a*rho;
  position(e,:) = 1/(6*a)*[sum((vx+vxS).*temp),sum((vy+vyS).*temp)];
end

volume0 = volume;

for p=1:pCount
  data       = getCPDIPolygonData(p,particles,mesh);
  nodeid{p}  = data.nodes;
  funcW{p}   = data.wf;
  gradW{p,1} = data.wg(:,1);
  gradW{p,2} = data.wg(:,2);
end

%% node quantities

nmass     = zeros(mesh.nodeCount,1);  % nodal mass vector
nmomentum = zeros(mesh.nodeCount,2);  % nodal momentum vector
niforce   = zeros(mesh.nodeCount,2);  % nodal internal force vector
neforce   = zeros(mesh.nodeCount,2);  % nodal external force vector
nforce    = zeros(mesh.nodeCount,2);  % nodal force vector

%% New: store shape functions/grads/nodes for all particles
% to speed up 

basis    = cell(pCount,1);
grad     = cell(pCount,1);

%% plot mesh, particles

%plot(node(fixNodes,1),node(fixNodes,2),'r*','markersize',14);

figure
hold on
MaxNVer = max(cellfun(@numel,Element));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Element,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Node,'FaceColor','w'); pause(1e-6)
%patch('Faces',ElemMat(trackedParId,:),'Vertices',Node,'FaceColor','r'); pause(1e-6)
plot_mesh(node,element,'Q4','k-',1.6); % background grid
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(position(:,1),position(:,2),'k.','markersize',10);
plot(Node(markedNode,1),Node(markedNode,2),'r.','markersize',25);
plot(node(fixNodes,1),node(fixNodes,2),'b*','markersize',10);
axis on % axis off

ta = 0;           % time
ka = 0;

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-16; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.2*(mesh.deltax/c);
time  = 3.;
t     = 0;
nsteps = floor(time/dtime);

istep = 1;

while ( t < time )
    disp(['time step ',num2str(t)])
    % reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    neforce(:)   = 0;
    % loop over particles
    for p=1:pCount
        sig    = stress(p,:);
        % particle mass and momentum to node        
        input.nodes=nodeid{p};
        input.wf   = funcW{p};
        input.wg   = [gradW{p,1} gradW{p,2}];
        input.Vp   = volume(p);
        
        shape = getCPDIPolygonBasis(p,input,particles,mesh);
        esctr = shape.node;
        % store basis function and gradients to speed up the code
        basis{p} = shape.phi;
        grad{p}  = shape.dphi;
        for i=1:length(esctr)
            id              = esctr(i);
            nmass(id)       = nmass(id)       + shape.phi(i)*mass(p);
            nmomentum(id,:) = nmomentum(id,:) + shape.phi(i)*mass(p)*velo(p,:);
            niforce(id,1)   = niforce(id,1) - volume(p)*(sig(1)*shape.dphi(i,1) + sig(3)*shape.dphi(i,2));
            niforce(id,2)   = niforce(id,2) - volume(p)*(sig(3)*shape.dphi(i,1) + sig(2)*shape.dphi(i,2));            
            neforce(id,:)   = neforce(id,:) + mass(p)*shape.phi(i)*bodyf;           
        end
    end
        
    % update nodal momenta   
    nforce    = niforce + neforce;
    nforce   (fixNodes,:)  = 0;
    nmomentum(fixNodes,:)  = 0;
    
    nmomentum = nmomentum + nforce*dtime;
        
    % update particle velocity and position and stresses
    
    % loop over particles
    for p=1:pCount
        Lp    = zeros(2,2);            
        esctr = nodeid{p};
        phi   = basis{p};
        dphi  = grad{p};
        for i=1:length(esctr)
            id = esctr(i);
            vI = [0 0];
            if nmass(id) > tol
                velo(p,:)  = velo(p,:) + dtime * phi(i)*nforce(id,:)/nmass(id);                
                vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
            end
            Lp = Lp + vI'*dphi(i,:);         % particle gradient velocity
        end
        
        F          = ([1 0;0 1] + Lp*dtime)*reshape(deform(p,:),2,2);
        deform(p,:)= reshape(F,1,4);
        volume(p)  = det(F)*volume0(p);
        J       = det(F);
        b       = F*F';
        sigma   = 1/J*( mu*(b-I) + lambda*log(J)*I );
        stress(p,:)  = [sigma(1,1) sigma(2,2)  sigma(1,2) ];
    end
    
    % update particle corners positions
    for c=1:size(particles.node,1)        
        xc    = particles.node(c,:);
        ec    = point2ElemIndex(xc,mesh);
        esctr = element(ec,:);
        for i=1:length(esctr)
            id      = esctr(i);
            x       = xc - node(id,:);
            [N,dNdx]= getMPM2D(x,mesh.deltax,mesh.deltay);
            if nmass(id) > tol                
                xc = xc + dtime*N*nmomentum(id,:)/nmass(id);                
            end            
        end
        particles.node(c,:) = xc;
    end    
    
    % update CPDI data     
    for p=1:pCount
      data       = getCPDIPolygonData(p,particles,mesh);
      nodeid{p}  = data.nodes;
      funcW{p}   = data.wf;
      gradW{p,1} = data.wg(:,1);
      gradW{p,2} = data.wg(:,2);
      volum{p}   = data.Vp; 
      volume(p)  = sum(data.Vp);
    end
    
    % VTK output    
    if (  mod(istep,interval) == 0 )
        vtkFile = sprintf('../../results/cpdi/%s%d',vtkFileName,istep);
        data.stress  = [stress zeros(pCount,1)];
        %data.pstrain = [];
        data.color    = color;
        VTKPolygonParticles(particles,vtkFile,data);
    end
        
    % advance to the next time step    
    t     = t + dtime;
    istep = istep + 1;
    
    % store time,velocty for plotting
    u  = particles.node(markedNode,2)-markedYCoord;
    ta = [ta;t];
    ka = [ka;u];        
end

%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

% pvdFile = fopen(strcat('../results/cpdi2/verticalBar/',vtkFileName,'.pvd'), 'wt');
% 
% fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
% fprintf(pvdFile,'<Collection>\n');
% 
% for i = 1:nsteps
%     if (  mod(i,interval) == 0 )
%         vtuFile = sprintf('%s%d%s',vtkFileName,i,'.vtu');
%         fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
%     end
% end
% 
% fprintf(pvdFile,'</Collection>\n');
% fprintf(pvdFile,'</VTKFile>\n');
% 
% fclose(pvdFile);

%%
Ux= zeros(size(mesh.node,1),1);
Uy= zeros(size(mesh.node,1),1);
sigmaXX = zeros(size(mesh.node,1),1);
sigmaYY = zeros(size(mesh.node,1),1);
sigmaXY = zeros(size(mesh.node,1),1);

VTKPostProcess(mesh.node,mesh.element,2,'Quad4',vtkFileName1,...
    [sigmaXX sigmaYY sigmaXY],[Ux Uy]);

%%

ss=load('../../example3D/cpdi/cpdiTet4Cantilever.mat');

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
plot(ss.ta(1:end),ss.ka(1:end),'r-','LineWidth',1.6);
xlabel('Time (s)')
ylabel('Displacement [m]')
legend('CPDI-T3','CPDI-T4')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
%axis([0 0.25 200 -1600])
box on 
%%
figure
hold on
plot_mesh(node,element,'Q4','k-',1.); % background grid
patch('Faces',ElemMat,'Vertices',particles.node,'FaceColor','w'); pause(1e-6)
%patch('Faces',ElemMat,'Vertices',Node,'FaceColor','b'); pause(1e-6)
%patch('Faces',ElemMat(trackedParId,:),'Vertices',Node,'FaceColor','r'); pause(1e-6)
%plot(coords3(1),coords3(2),'r.','markersize',15);

disp([num2str(toc),'   DONE '])
