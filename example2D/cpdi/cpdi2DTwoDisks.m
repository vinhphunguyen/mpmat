% This file implements the Material Point Method with CPDI interpolation
% described in the article
%
% A. Sadeghirad, R. M. Brannon, and J. Burghardt. A convected particle domain
% interpolation technique to extend applicability of the material point
% method for problems involving massive deformations. IJNME, 86(12):1435--1456, 2011.
%
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
%
% Two elastic disks come into contact.
% Example taken from Sulsky et al. 1994 paper.
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

E   = 1000;        % Young's modulus
nu  = 0.3;         % Poisson ratio
rho = 1000;        % density
kappa = 3-4*nu;    % Kolosov constant
mu    = E/2/(1+nu);% shear modulus

v   = 0.1;

interval     = 100;
vtkFileName  = 'cpdi2DTwoDisks';

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell=1;
l     = 1;
numx2 = 20;      % number of elements along X direction
numy2 = 20;      % number of elements along Y direction
[mesh]= buildGrid2D(l,l,numx2,numx2, ghostCell);
element= mesh.element;
node   = mesh.node;

%% particles
noParticleX = 2;
noParticleY = 2;
% noParticlexnoParticle Gaussian quadrature
%[W,Q]=quadrature( noParticleX, 'GAUSS', 2 );
lpx    = mesh.deltax/noParticleX;
lpy    = mesh.deltay/noParticleY;
w      = lpx*lpy;
% body1

center = [0.2+mesh.deltax*ghostCell 0.2+mesh.deltay*ghostCell];
radius = 0.2;

vo = []; ma   = []; co = [];

for e=1:mesh.elemCount                          % start of element loop
    sctr = element(e,:);                        %  element scatter vector
    pts  = node(sctr,:);
    for j=1:noParticleY
        for i=1:noParticleX
            x(1)   = pts(1,1) + lpx*0.5 + (i-1)*lpx;
            x(2)   = pts(1,2) + lpy*0.5 + (j-1)*lpy;
            r      = norm(x-center);
            if ( r-radius < 0 )
                vo  = [vo;w];
                ma  = [ma; w*rho];
                co  = [co;x];
            end
        end
    end
end

% body 2
center = [0.8+mesh.deltax*ghostCell 0.8+mesh.deltay*ghostCell];
radius = 0.2;

for e=1:mesh.elemCount                     % start of element loop
    sctr = element(e,:);                   %  element scatter vector
    pts  = node(sctr,:);
    for j=1:noParticleY
        for i=1:noParticleX
            x(1)   = pts(1,1) + lpx*0.5 + (i-1)*lpx;
            x(2)   = pts(1,2) + lpy*0.5 + (j-1)*lpy;
            r      = norm(x-center);
            if ( r-radius < 0 )
                vo  = [vo;w];
                ma  = [ma; w*rho];
                co  = [co;x];
            end
        end
    end
end

pCount = length(vo);                       % # of particles
volume  = vo;
mass    = ma;
coord   = co;
volume0 = volume;
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stressx
strain  = zeros(pCount,3);                % strain
velo    = [ones(pCount/2,2)*v;-ones(pCount/2,2)*v];% velocity
dvec1   = zeros(pCount,2);                % domain vector 1, r1
dvec2   = zeros(pCount,2);                % domain vector 2, r2

for p=1:pCount
    dvec1(p,1) = 0.5*lpx;
    dvec1(p,2) = 0;
    dvec2(p,1) = 0;
    dvec2(p,2) = 0.5*lpy;
end

dvec10 = dvec1;
dvec20 = dvec2;

%% node quantities

nmass     = zeros(mesh.nodeCount,1);  % nodal mass vector
nmomentum = zeros(mesh.nodeCount,2);  % nodal momentum vector
niforce   = zeros(mesh.nodeCount,2);  % nodal internal force vector
neforce   = zeros(mesh.nodeCount,2);  % nodal external force vector

%% plot mesh, particles

hold on
%plot_mesh(node1,element1,elemType,'r-',0.2);
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(coord(:,1),coord(:,2),'k.','markersize',10);
axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.001;
time  = 3.5;
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
        xp     = coord(p,:);
        r1p    = dvec1(p,:);
        r2p    = dvec2(p,:);
        % particle mass and momentum to node
        data  = cpdi2D(xp,r1p,r2p,mesh);
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
        xp    = coord(p,:);
        r1p   = dvec1(p,:);
        r2p   = dvec2(p,:);
        data  = cpdi2D(xp,r1p,r2p,mesh);
        esctr = data.node;        
        for i=1:length(esctr)
            id = esctr(i);
            vI = [0 0];
            if nmass(id) > tol
                velo(p,:)  = velo(p,:) + dtime * data.phi(i)*niforce(id,:)  /nmass(id);
                coord(p,:) = coord(p,:)+ dtime * data.phi(i)*nmomentum(id,:)/nmass(id);
                vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
            end
            Lp = Lp + vI'*data.dphi(i,:);         % particle gradient velocity
        end
        
        F          = ([1 0;0 1] + Lp*dtime)*reshape(deform(p,:),2,2);
        deform(p,:)= reshape(F,1,4);
        dvec1(p,:) = dvec10(p,:)*F';
        dvec2(p,:) = dvec20(p,:)*F';
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
        
    % store time,velocty for plotting
    
    ta = [ta;t];
    ka = [ka;k];
    sa = [sa;u];
    
    % VTK output
    
    if (  mod(istep,interval) == 0 )        
        vtkFile = sprintf('../results/cpdi/%s%d',vtkFileName,istep);
        data.stress  = [stress zeros(pCount,1)];
        data.pstrain = [];
        data.velo    = velo;
        VTKParticles(coord,vtkFile,data);
    end
    
    
    % advance to the next time step
    
    t     = t + dtime;
    istep = istep + 1;
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

pvdFile = fopen(strcat('../results/',vtkFileName,'.pvd'), 'wt');

fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
fprintf(pvdFile,'<Collection>\n');

for i = 1:nsteps
    if (  mod(i,interval) == 0 )
        vtuFile = sprintf('%s%d%s',vtkFileName,i,'.vtp');
        fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
    end
end

fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
plot(ta(1:end),sa(1:end),'r--','LineWidth',2);
plot(ta(1:end),ka(1:end)+sa(1:end),'g-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy')
legend('kinetic','strain','total')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
axis([0 3 0 3])

ss=load('gimp-2disks-energies.mat');

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'g-','LineWidth',2.1);
plot(ss.ta(1:end),ss.ka(1:end),'r-','LineWidth',2.1);
plot(ta(1:end),sa(1:end),'g--','LineWidth',2.1);
plot(ss.ta(1:end),ss.sa(1:end),'r--','LineWidth',2.1);
xlabel('Time')
ylabel('Energy')
legend('cpdi','gimp')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])

disp([num2str(toc),'   DONE '])
