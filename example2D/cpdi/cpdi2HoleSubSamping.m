% This file implements the Material Point Method with CPDI-Q4 interpolation
% described in the article
%
% A. Sadeghirad, R. M. Brannon, and J.E. Guilkey. Second-order convected
% particle domain in- terpolation (CPDI2) with enrichment for weak discontinuities
% at material interfaces. IJNME, 95(11):928-952, 2013.
%
% This example demonstrates that holes can appear in CPDI. And solution
% using the sub-sampling technique.
%
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia
% 14 August 2015.

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

E   = 1000;        % Young's modulus
nu  = 0.3;         % Poisson ratio
rho = 1000;        % density
kappa = 3-4*nu;    % Kolosov constant
mu    = E/2/(1+nu);% shear modulus

v   = 0.5;

interval     = 1;
vtkFileName  = 'cpdi2Hole';
vtkFileName1 = '../results/cpdi2/cpdi2HoleGrid';

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell=0;
l     = 7;
numx2 = 7;      % number of elements along X direction
numy2 = 7;      % number of elements along Y direction
[mesh]= buildGrid2D(l,l,numx2,numx2, ghostCell);
element= mesh.element;
node   = mesh.node;

%%   particle distribution from a mesh
%

% this one causes hole
node1 = [1.5 2.5;
         5.5 2.5;
         5.5 4.5;
         1.5 4.5;
         3.5 2.5; % extra nodes for sub-sampling
         3.5 4.5; % extra nodes for sub-sampling
         %% particle 2
         3.1 6.1;
         3.9 6.1;
         3.9 6.9;
         3.1 6.9];

% this one does not
% node1 = [2.5 2.5;
%          4.5 2.5;
%          4.5 4.5;
%          2.5 4.5;
%          3.1 6.1;
%          3.9 6.1;
%          3.9 6.9;
%          3.1 6.9];       

element1    = cell(2,1);
element1{1} = [1 5 6 4 5 2 3 6]; 
element1{2} = [7 8 9 10];      
          
numelem  = size(element1,1);
numnode  = size(node1,1); 
elemType = 'Q4';

% store the particle mesh into a structure for convenience
particles.node     = node1;
particles.elem     = element1;
particles.elemType = elemType;

pCount  = numelem;                       % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
coords  = zeros(pCount,2);
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stressx
strain  = zeros(pCount,3);                % strain
velo    = zeros(pCount,2);
% CPDI data
nodeid  = cell(pCount,1);                 % nodes affect particle 'p'
funcW   = cell(pCount,1);                 % function weights of 'p'
gradW   = cell(pCount,2);                 % gradient weights of 'p', 
color   = zeros(pCount,1);                % 1 column for grad_x 
gradW   = cell(pCount,2);
volum   = cell(pCount,1);                 % store particle domain volume

color(1) = 1;
color(2) = 2;

% particle mass, volume and initial velocity
for e = 1:numelem
    coord = node1(element1{e},:);
    a     = 0.5*( coord(1,1)*coord(2,2)  - coord(2,1)*coord(1,2) ...
                + coord(2,1)*coord(3,2)  - coord(3,1)*coord(2,2) ...
                + coord(3,1)*coord(4,2)  - coord(4,1)*coord(3,2) ...
                + coord(4,1)*coord(1,2)  - coord(1,1)*coord(4,2) );
    volume(e)   = a;
    mass(e)     = a*rho;
    coords(e,:) = mean(coord); % center of each element=particle
end

volume(1) = volume(1);
%mass(1)   = mass(1)*2;
volum{1} = [volume(1) volume(1)];
volum{2} = volume(2);

velo(2,2) = -v;
volume0 = volume;

for p=1:pCount
  data       = getCPDIQuadDataGeneral(p,particles,mesh);
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

%% plot mesh, particles

% hold on
% plot_mesh(particles.node,particles.elem,elemType,'b-',0.2);
% plot_mesh(particles.node,particles.elem(1,:),elemType,'r--',0.2);
% plot_mesh(node,element,'Q4','k-',1.);
% %plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
% plot(coords(:,1),coords(:,2),'k.','markersize',10);
% axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.05*mesh.deltax/c;
time  = 10.5;
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
    disp(['time step ',num2str(t)])
    % reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    % loop over particles
    for p=1:pCount
        sig    = stress(p,:);
        % particle mass and momentum to node
        %data  = cpdi22D(p,particles,mesh);
        input.nodes=nodeid{p};
        input.wf   = funcW{p};
        input.wg   = [gradW{p,1} gradW{p,2}];
        input.Vp   = volume(p);
        
        data  = getCPDIQuadBasisGeneral(p,input,particles,mesh);
        esctr = data.node;
        for i=1:length(esctr)
            id              = esctr(i);
            nmass(id)       = nmass(id)       + data.phi(i)*mass(p);
            nmomentum(id,:) = nmomentum(id,:) + data.phi(i)*mass(p)*velo(p,:);
            niforce(id,1)   = niforce(id,1) - volume(p)*(sig(1)*data.dphi(i,1) + sig(3)*data.dphi(i,2));
            niforce(id,2)   = niforce(id,2) - volume(p)*(sig(3)*data.dphi(i,1) + sig(2)*data.dphi(i,2));
        end
    end
    
    % debug
    % update nodal momenta
    
    nmomentum = nmomentum + niforce*dtime;
    % no bounda
    % update particle velocity and position and stresses
    k = 0;
    u = 0;
    % loop over particles
    for p=1:pCount
        Lp    = zeros(2,2);        
        %data  = cpdi22D(p,particles,mesh); % old implementation
        input.nodes=nodeid{p};
        input.wf   = funcW{p};
        input.wg   = [gradW{p,1} gradW{p,2}];
        input.Vp   = volume(p);
        
        data  = getCPDIQuadBasisGeneral(p,input,particles,mesh);
        esctr = data.node;
        for i=1:length(esctr)
            id = esctr(i);
            vI = [0 0];
            if nmass(id) > tol
                velo(p,:)  = velo(p,:) + dtime * data.phi(i)*niforce(id,:)  /nmass(id);
                %coords(p,:)= coords(p,:)+ dtime * data.phi(i)*nmomentum(id,:)/nmass(id);
                vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
            end
            Lp = Lp + vI'*data.dphi(i,:);         % particle gradient velocity
        end
        
        F          = ([1 0;0 1] + Lp*dtime)*reshape(deform(p,:),2,2);
        deform(p,:)= reshape(F,1,4);
        volume(p)  = det(F)*volume0(p);
        dEps       = dtime * 0.5 * (Lp+Lp');
        dsigma     = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
        stress(p,:)= stress(p,:) + dsigma';
        strain(p,:)= strain(p,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
        
        k = k + 0.5*(velo(p,1)^2+velo(p,2)^2)*mass(p);
        %             u = u + 0.25/mu*( 0.25*(kappa+1)*(s(pid,1)^2+s(pid,2)^2) ...
        %                   - 2*(s(pid,1)*s(pid,2)-s(pid,3)^2))*Vp(pid);
        u = u + 0.5*volume(p)*stress(p,:)*strain(p,:)';
    end
    
    % update particle corners position
    for c=1:numnode        
        xc    = particles.node(c,:);
        ec    = point2ElemIndex(xc,mesh);
        esctr = element(ec,:);
        for i=1:length(esctr)
            id = esctr(i);
            x     = xc - node(id,:);
            [N,dNdx]=getMPM2D(x,mesh.deltax,mesh.deltay);
            if nmass(id) > tol                
                xc = xc + dtime*N*nmomentum(id,:)/nmass(id);                
            end            
        end
        particles.node(c,:) = xc;
    end
    
    % update CPDI data 
    
    for p=1:pCount
      data = getCPDIQuadDataGeneral(p,particles,mesh);
      nodeid{p}  = data.nodes;
      funcW{p}   = data.wf;
      gradW{p,1} = data.wg(:,1);
      gradW{p,2} = data.wg(:,2);
      volum{p}   = data.Vp; 
      volume(p)=sum(data.Vp);
    end
    
    % store time,velocty for plotting
    
    ta = [ta;t];
    ka = [ka;k];
    sa = [sa;u];
    
    % VTK output
    
    if (  mod(istep,interval) == 0 )
        vtkFile = sprintf('../results/cpdi2/%s%d',vtkFileName,istep);
        data.stress  = [stress zeros(pCount,1)];
        data.pstrain = [];
        data.velo    = velo;
        data.color   = color;
        VTKParticlesCPDI(particles,vtkFile,data);
    end
    
    
    % advance to the next time step
    
    t     = t + dtime;
    istep = istep + 1;
end



%% post processing

% disp([num2str(toc),'   POST-PROCESSING '])
% 
% pvdFile = fopen(strcat('../results/cpdi2/',vtkFileName,'.pvd'), 'wt');
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
% Ux= zeros(size(mesh.node,1),1);
% Uy= zeros(size(mesh.node,1),1);
% sigmaXX = zeros(size(mesh.node,1),1);
% sigmaYY = zeros(size(mesh.node,1),1);
% sigmaXY = zeros(size(mesh.node,1),1);
% 
% VTKPostProcess(mesh.node,mesh.element,2,'Quad4',vtkFileName1,...
%     [sigmaXX sigmaYY sigmaXY],[Ux Uy]);

%%

ss=load('cpdi-hole-normal-energies.mat');

figure
subplot(2,1,1)
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'black-','LineWidth',2.1);
plot(ss.ta(1:end),ss.ka(1:end),'r-','LineWidth',2.1);
xlabel('Time')
ylabel('Kinetic energy')
legend('sub-sampling','standard')

% 
subplot(2,1,2)
set(gca,'FontSize',14)
hold on
plot(ta(1:end),sa(1:end),'black-','LineWidth',2.1);
plot(ss.ta(1:end),ss.sa(1:end),'r-','LineWidth',2.1);
xlabel('Time')
ylabel('Strain energy')
legend('sub-sampling','standard')

% figure
% hold on
% plot_mesh(particles.node,particles.elem,elemType,'r-',0.2);
% plot_mesh(node,element,'Q4','k-',1.);
% %plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
% plot(coords(:,1),coords(:,2),'k.','markersize',10);
% axis off

disp([num2str(toc),'   DONE '])
