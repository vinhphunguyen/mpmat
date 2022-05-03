% This file implements the Material Point Method of Sulsky 1994.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% Particles are generated using an unstructured FE mesh (gmsh).
%
% Shape functions: using standard FE shape function routine by converting
% particle position to natural coordinates.
%
% Implicit time integration described in
%
% D. Sulsky and A. Kaul. Implicit dynamics in the material-point method.
% Computer Methods in Applied Mechanics and Engineering, 193(12-14):1137?1170, 2004.
%
%
% Two elastic disks come into contact.
% Example taken from Sulsky et al. 1994 paper.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% June 2014.

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
v     = 0.1;     % initial particle velocity
bulk  = E/3/(1-2*nu);

interval     = 100;% time interval for saving vtp files.
vtkFileName  = 'mpm2DTwoDisks';

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])

%%   particle distribution from a mesh
%

meshFile = 'disks.msh';
mesh     = load_gmsh (meshFile);

elemType = 'T3';
numnode  = mesh.nbNod;
numelem  = mesh.nbTriangles;
node1    = mesh.POS(:,1:2);
element1 = mesh.TRIANGLES(1:numelem,1:3);

% check if Jacobian is negative

element1  = tricheck(node1,element1,1);

% one particle per triangle

pCount = numelem;

mass    = ones(pCount,1);                 % mass
volume  = ones(pCount,1);                 % volume
deform  = ones(pCount,4);                 % gradient deformation
stress  = zeros(pCount,3);                % stress
strain  = zeros(pCount,3);                % strain
velo    = zeros(pCount,2);                % velocity
coords  = zeros(pCount,2);                % position

% initialise particle position, mass ,volume, velocity

for e = 1:numelem
    coord = node1(element1(e,:),:);
    a     = det([coord,[1;1;1]])/2;
    volume(e)   = a;
    mass(e)     = a*rho;
    coords(e,:) = mean(coord); % center of each element=particle
    
    if coords(e,1) < 0.5
        velo(e,:) = [v v];
    else
        velo(e,:) = [-v -v];
    end
    
    deform(e,:) = [1 0 0 1];
end
volume0 = volume;

%% Computational grid

l = 1;

numx2 = 20;      % number of elements along X direction
numy2 = 20;      % number of elements along Y direction

deltax = l/numx2;
deltay = l/numy2;

nnx=numx2+1;
nny=numy2+1;
node=square_node_array([0 0],[l 0],[l l],[0 l],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element =make_elem(node_pattern,numx2,numy2,inc_u,inc_v);

elemCount = size(element,1);
nodeCount = nnx*nny;

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

pElems  = ones(pCount,1);
mpoints = cell(elemCount,1);

for p=1:pCount
    x = coords(p,1);
    y = coords(p,2);
    e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
    pElems(p) = e;
end

for e=1:elemCount
    id  = find(pElems==e);
    mpoints{e}=id;
end

%% node quantities

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,2);  % nodal momentum vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector (no need for this exam)
K         = sparse(nodeCount*2,nodeCount*2); % stiffness matrix
M         = sparse(nodeCount*2,nodeCount*2); % mass matrix
activeNodes = unique(element(unique(pElems),:));
activeDofs(1:2:2*length(activeNodes)) = activeNodes*2-1;
activeDofs(2:2:2*length(activeNodes)) = activeNodes*2;
nvelo     = zeros(2*nodeCount,1);
nvelo0    = zeros(2*nodeCount,1);

%% plot mesh, particles

figure(1)
hold on
plot_mesh(node1,element1,elemType,'r-',0.2);
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(coords(:,1),coords(:,2),'k.','markersize',10);
axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.01;
time  = dtime;
t     = 0;

dtime2 = dtime*dtime;
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
    K(:,:)       = 0;
    % loop over computational cells or elements
    for e=1:elemCount
        esctr = element(e,:);      % element connectivity
        enode = node(esctr,:);     % element node coords
        mpts  = mpoints{e};        % particles inside element e
        for p=1:length(mpts)       % loop over particles
            pid  = mpts(p);
            % particle mass and momentum to node
            sigma  = stress(pid,:);
            vol    = volume(pid);
            m      = mass(pid);
            vel    = velo(pid,:);
            xpa    = coords(pid,:);
            
            pt(1)= (2*xpa(1)-(enode(1,1)+enode(2,1)))/deltax;
            pt(2)= (2*xpa(2)-(enode(2,2)+enode(3,2)))/deltay;
            
            [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
            J0       = enode'*dNdxi;             % element Jacobian matrix
            invJ0    = inv(J0);
            dNdx     = dNdxi*invJ0;
            
            for i=1:length(esctr)
                id    = esctr(i);    
                BI      = [dNdx(i,1) 0;0 dNdx(i,2);dNdx(i,2) dNdx(i,1)];
                nmass(id)       = nmass(id)       + N(i)*m;
                nmomentum(id,:) = nmomentum(id,:) + N(i)*m*vel;
                niforce(id,:)   = niforce(id,:) - vol*sigma*BI;
                for j=1:length(esctr)
                    jd    = esctr(j);
                    BJ      = [dNdx(j,1) 0;0 dNdx(j,2);dNdx(j,2) dNdx(j,1)];
                    K([2*id-1 2*id],[2*jd-1 2*jd]) = K([2*id-1 2*id],[2*jd-1 2*jd]) + ...
                        vol*BI'*C*BJ;
                end
            end
        end
    end
    
    nvelo0(2*activeNodes-1) = nmomentum(activeNodes,1)./nmass(activeNodes);
    nvelo0(2*activeNodes  ) = nmomentum(activeNodes,2)./nmass(activeNodes);
    
    for i=1:nodeCount
        M(2*i-1,2*i-1) = nmass(i);
        M(2*i,2*i)     = nmass(i);
    end
    
    % update nodal velocity by solving Ax=b
    
    nmomentum = nmomentum + niforce*dtime;
    A         = (M+dtime2*K);
    f         = reshape(nmomentum',2*nodeCount,1);
    
    % apply boundary conditions
    freeNodes= setdiff(1:nodeCount,activeNodes)';
    udofs    = [2*freeNodes-1];
    vdofs    = [2*freeNodes];
    uFixed   = [zeros(length(freeNodes),1)];
    vFixed   = [zeros(length(vdofs),1)];
    [A,f]    = applyDirichletBCs(A,f,udofs,vdofs,uFixed',vFixed');
    % solve the system
    nvelo = A\f;
    
    %nvelo(activeDofs) = A(activeDofs,activeDofs)\f(activeDofs);
         
    %nvelo     = A\f;
    
    acce      = (nvelo-nvelo0);
    
    % update particle velocity and position and stresses
    k = 0;
    u = 0;
    for e=1:elemCount
        esctr = element(e,:);
        enode = node(esctr,:);
        mpts  = mpoints{e};        
        for p=1:length(mpts) % loop over particles
            pid  = mpts(p);
            xpa  = coords(pid,:);
            Lp   = zeros(2,2);
            for i=1:length(esctr)
                id = esctr(i);
                vI = nvelo([2*id-1;2*id]);
                aI = acce([2*id-1;2*id]);
                x       = xpa - node(id,:);
                [N,dNdx]= getMPM2D(x,deltax,deltay);                                
                velo(pid,:)  = velo(pid,:)   +  N*aI';
                coords(pid,:)= coords(pid,:) + dtime * N*vI';                
                Lp = Lp + vI*dNdx;
            end
            
            F              = ([1 0;0 1] + Lp*dtime)*reshape(deform(pid,:),2,2);
            deform(pid,:)  = reshape(F,1,4);
            volume(pid)    = det(F)*volume0(pid);
            dEps           = dtime * 0.5 * (Lp+Lp');
            dsigma         = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
            stress(pid,:)  = stress(pid,:) + dsigma';
            strain(pid,:)  = strain(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
            
            k = k + 0.5*(velo(pid,1)^2+velo(pid,2)^2)*mass(pid);   
            u = u + 0.5*volume(pid)*stress(pid,:)*strain(pid,:)';
        end
    end
    
    % store time,velocty for plotting
%     
%     pos{istep} = coords;
%     vel{istep} = velo;
    
    % update the element particle list
    
    for p=1:pCount
        x = coords(p,1);
        y = coords(p,2);
        e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
        pElems(p) = e;
    end
    
    for e=1:elemCount
        id  = find(pElems==e);
        mpoints{e}=id;
    end
    
    activeNodes         = unique(element(unique(pElems),:));
    activeDofs=[];
    activeDofs(1:2:2*length(activeNodes)) = activeNodes*2-1;
    activeDofs(2:2:2*length(activeNodes)) = activeNodes*2;

    % store time,velocty for plotting
    
    ta = [ta;t];
    ka = [ka;k];
    sa = [sa;u];
    
    % VTK output
    
    %     if (  mod(istep,interval) == 0 )
    %         xp = pos{istep};
    %         vtkFile = sprintf('../results/%s%d',vtkFileName,istep);
    %         VTKParticles(xp,vtkFile,s);
    %     end
    
    
    % advance to the next time step
    
    t     = t + dtime;
    istep = istep + 1;
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

ss=load('mpmTwoDisksImplicit5.mat');
% 
ka = ka*1000;
sa = sa*1000;

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'black-','LineWidth',1.6);
plot(ta(1:end),sa(1:end),'r-','LineWidth',2);
plot(ta(1:end),ka(1:end)+sa(1:end),'g-','LineWidth',2.1);
plot(ss.ta(1:end),ss.ka(1:end),'black-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy')
legend('kinetic','strain','total')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
%axis([0 3 0 3])

%savefile = 'mpm-implicit-2disks-dt002.mat';
%save(savefile,'ta','ka','sa');

disp([num2str(toc),'   DONE '])

%% plot curves for different time steps

ss1 = load('mpm-implicit-2disks-dt001.mat');
ss2 = load('mpm-implicit-2disks-dt002.mat');

figure
set(gca,'FontSize',14)
hold on
plot(ss1.ta(1:end),ss1.ka(1:end),'black-','LineWidth',2.1);
plot(ss1.ta(1:end),ss1.sa(1:end),'red-','LineWidth',2.1);
plot(ss1.ta(1:end),ss1.ka(1:end)+ss1.sa(1:end),'g-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy')
legend('kinetic','strain','total')

exportfig(gcf,'implicit-sulsky-2disks-dt01.eps',opts)


