% This file implements the Material Point Method with CPDI2 interpolation
% described in the article (but with triangular particle domains)
%
% Coined as CPDI-T3 in my book.
%
% A. Sadeghirad, R. M. Brannon, and J.E. Guilkey. Second-order convected
% particle domain in- terpolation (CPDI2) with enrichment for weak discontinuities
% at material interfaces. IJNME, 95(11):928-952, 2013.
%
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% Particles taken from a T3 mesh (Gmsh).
%
% Vertical bar with extreme deformation.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% May 2014, Saigon, Vietnam.

%%

addpath ../fem_util/
addpath ../fem-functions/
addpath ../post-processing/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
%

E   = 1;           % Young's modulus
nu  = 0.3;         % Poisson ratio
rho = 1050e-12;    % density
K   = E/3/(1-2*nu);    % bulk modulus
mu    = E/2/(1+nu);% shear modulus
lambda = K - 2/3*mu;

g     = 1000e3; % gravity
bodyf = [0 -g];

I  = [1 0;0 1];

interval     = 1;
vtkFileName  = 'bar';
vtkFileName1 = '../results/cpdi3/verticalBar/barGrid';

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])

%%   particle distribution
%

l     = 1000;
numx2 = 7;      % number of elements along X direction
numy2 = 7;      % number of elements along Y direction

[node,elemLst]=make_cross_mesh([0 0],[l l],numx2,numy2);

node(:,1) = node(:,1) +   l/2;
node(:,2) = node(:,2) + 4*l/2;

% store the particle mesh into a structure for convenience
elemType           = 'T3';
particles.node     = node;
particles.elem     = elemLst{:,5};
particles.elemType = elemType;

pCount  = size(particles.elem,1);         % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
coords  = zeros(pCount,2);
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stressx
strain  = zeros(pCount,3);                % strain
velo    = zeros(pCount,2);


for e = 1:pCount
    coord = particles.node(particles.elem(e,:),:);
    a     = 0.5*( (coord(1,2)-coord(3,2))*(coord(1,1)-coord(2,1))...
                 +(coord(1,2)-coord(2,2))*(coord(3,1)-coord(1,1)) );
    volume(e)  = a;
    mass(e)    = a*rho;
    coords(e,:) = mean(coord); % center of each element=particle  
end

volume0 = volume;

%% Computational grid

ghostCell=0;
lx     = (l/2)*4;
ly     = (l/2)*7;
numx2 = 4;      % number of elements along X direction
numy2 = 7;      % number of elements along Y direction
[mesh]= buildGrid2D(lx,ly,numx2,numy2, ghostCell);
element= mesh.element;
node   = mesh.node;

% find boundary nodes

fixNodes=find(abs(node(:,2)-6*l/2)<1e-10);

%% node quantities

nmass     = zeros(mesh.nodeCount,1);  % nodal mass vector
nmomentum = zeros(mesh.nodeCount,2);  % nodal momentum vector
niforce   = zeros(mesh.nodeCount,2);  % nodal internal force vector
neforce   = zeros(mesh.nodeCount,2);  % nodal external force vector
nforce    = zeros(mesh.nodeCount,2);  % nodal force vector

%% plot mesh, particles

hold on
plot_mesh(particles.node,particles.elem,elemType,'black-',1.2);
plot_mesh(node,element,'Q4','r-',1.2);
plot(coords(:,1),coords(:,2),'b*','markersize',10);
%plot(node(fixNodes,1),node(fixNodes,2),'r*','markersize',14);
axis off

ta = 0;           % time
ka = 0;

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-16; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.1*(mesh.deltax/c);
time  = 0.25;
t     = 0;
nsteps = floor(time/dtime);

pNodeOut = cell(nsteps,1);
pElemOut = cell(nsteps,1);

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
        data  = cpdi2Triangle(p,particles,mesh);
        esctr = data.node;
        for i=1:length(esctr)
            id              = esctr(i);
            nmass(id)       = nmass(id)       + data.phi(i)*mass(p);
            nmomentum(id,:) = nmomentum(id,:) + data.phi(i)*mass(p)*velo(p,:);
            niforce(id,1)   = niforce(id,1) - volume(p)*(sig(1)*data.dphi(i,1) + sig(3)*data.dphi(i,2));
            niforce(id,2)   = niforce(id,2) - volume(p)*(sig(3)*data.dphi(i,1) + sig(2)*data.dphi(i,2));            
            neforce(id,:)   = neforce(id,:) + mass(p)*data.phi(i)*bodyf;           
        end
    end
        
    % update nodal momenta   
    nforce    = niforce + neforce;
    nforce   (fixNodes,2)  = 0;
    nmomentum(fixNodes,2)  = 0;
    
    nmomentum = nmomentum + nforce*dtime;
        
    % update particle velocity and position and stresses
    
    % loop over particles
    for p=1:pCount
        Lp    = zeros(2,2);        
        data  = cpdi2Triangle(p,particles,mesh);
        esctr = data.node;
        for i=1:length(esctr)
            id = esctr(i);
            vI = [0 0];
            if nmass(id) > tol
                velo(p,:)  = velo(p,:) + dtime * data.phi(i)*nforce(id,:)/nmass(id);                
                vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
            end
            Lp = Lp + vI'*data.dphi(i,:);         % particle gradient velocity
        end
        
        F          = ([1 0;0 1] + Lp*dtime)*reshape(deform(p,:),2,2);
        deform(p,:)= reshape(F,1,4);
        volume(p)  = det(F)*volume0(p);
        J       = det(F);
        b       = F*F';
        sigma   = 1/J*( mu*(b-I) + lambda*log(J)*I );
        stress(p,:)  = [sigma(1,1) sigma(2,2)  sigma(1,2) ];
    end
    % update particle corner position
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
    
    pNodeOut{istep} = particles.node;
      
    % VTK output
    
    if (  mod(istep,interval) == 0 )
        vtkFile = sprintf('../results/cpdi3/verticalBar/%s%d',vtkFileName,istep);
        data.stress  = [stress zeros(pCount,1)];
        data.pstrain = [];
        data.velo    = velo;
        VTKParticlesCPDI(particles,vtkFile,data);
    end
    
    % advance to the next time step
    
    t     = t + dtime;
    istep = istep + 1;
    
      % store time,velocty for plotting
    coord = particles.node(particles.elem(3,:),:);
    coords3 = mean(coord); % center of each element=particle
    u  = coords3(:,2)-coords(3,2);
    ta = [ta;t];
    ka = [ka;u];        
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

% pvdFile = fopen(strcat('../results/cpdi3/verticalBar/',vtkFileName,'.pvd'), 'wt');
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

ss=load('cpdi2VerticalBar.mat');

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
plot(ss.ta(1:end),ss.ka(1:end),'r-','LineWidth',1.6);
xlabel('Time')
ylabel('Displacement')
legend('CPDI2Tri','CPDI2Quad')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
axis([0 0.25 0 -1800])

step = 60;

pNodes = pNodeOut{step};

for e = 1:pCount
    coord = pNodes(particles.elem(e,:),:);
    coords(e,:) = mean(coord); % center of each element=particle
end

figure
hold on
plot_mesh(pNodes,particles.elem,elemType,'r-',1.05);
plot_mesh(node,element,'Q4','k-',1.2)
plot(coords(:,1),coords(:,2),'k.','markersize',10);
axis off

disp([num2str(toc),'   DONE '])
