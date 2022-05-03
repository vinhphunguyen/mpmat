% This file implements the Material Point Method of Sulsky 1994.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
%
%
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia
% 25 August 2015.

%%

addpath ../../grid/
addpath ../../basis/
addpath ../../particleGen/
addpath ../../constitutiveModels/
addpath ../../util/
addpath ../../postProcessing/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties 
% 

E   = 0.073e9;     % Young's modulus
nu  = 0.4;         % Poisson ratio
rho = 1.01e3;        % density
kappa = 3-4*nu;    % Kolosov constant
mu    = E/2/(1+nu);% shear modulus

v0    = 20;      % initial particle velocity [m/s]

interval     = 20;% time interval for saving vtp files.
vtkFileName  = 'mpm2DDiskWallNoSlip';

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid
lx = 12e-2;
ly = 12e-2;
ghostCell=0;
numx2 = 40;      % number of elements along X direction
numy2 = 40;      % number of elements along Y direction
[mesh]= buildGrid2D(lx,ly,numx2,numy2, ghostCell);
%mesh.node = mesh.node - [ones(mesh.nodeCount,1)*mesh.deltax ones(mesh.nodeCount,1)*mesh.deltay];

element   = mesh.element;
node      = mesh.node;
elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
deltax    = mesh.deltax;
deltay    = mesh.deltay;


% find boundary nodes

eps=1e-12;
left   = mesh.lNodes;
right  = mesh.rNodes;

%%   particle distribution 

circle1.center = (1/100)*[6 6];
circle2.center = circle1.center;
circle1.radius = 3e-2; % m
circle2.radius = 4e-2; % m

ppc = [2 2];

[res] = generateMPDiffCircles(circle1,circle2,ppc,mesh);
pCount = size(res.position,1);
body1.volume = res.volume;
body1.volume0 = res.volume;
body1.mass   = res.volume*rho;
body1.coord  = res.position;
body1.deform = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress = zeros(pCount,3);                % stress
body1.strain = zeros(pCount,3);                % strain
body1.velo   = zeros(pCount,2);               % velocity
body1.velo(:,1) = v0;

bodies(1) = body1;
bodyCount = length(bodies);

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles


for ib=1:length(bodies)
    body      = bodies(ib);
    elems     = ones(length(body.volume),1);
    
    for ip=1:length(body.volume)
        x = body.coord(ip,1); y = body.coord(ip,2);
        e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
        elems(ip) = e;
    end
    
    bodies(ib).elements = unique(elems);
    bodies(ib).nodes    = unique(element(bodies(ib).elements,:));
    
    mpoints = cell(elemCount,1);
    for ie=1:elemCount
        id  = find(elems==ie);
        mpoints{ie}=id;
    end    
    
    bodies(ib).mpoints  = mpoints;
end


%% plot mesh, particles

figure
hold on
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(bodies{1}.elements,:),'Q4','cy-',2.1);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','r-',2.1);
plot(bodies(1).coord(:,1),bodies(1).coord(:,2),'k.','markersize',10);
plot(node(left,1),node(left,2),'k.','markersize',10);
plot(node(right,1),node(right,2),'k.','markersize',10);
axis off

%% node quantities

nmassS    = zeros(nodeCount,1);  % nodal mass vector of the system
nmomentumS= zeros(nodeCount,2);  % nodal momentum vector of the system
niforceS  = zeros(nodeCount,2);  % nodal internal force of the system
nmass     = zeros(nodeCount,1);  % nodal mass vector of each body
nmomentum = zeros(nodeCount,2);  % nodal momentum vector of each body
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector
nvelo     = zeros(nodeCount,2);  % nodal velocities (body1,body2,center of mass)
nvelo0    = nvelo;


%%
disp([num2str(toc),'   SOLVING '])

tol   = 1e-12; % mass tolerance

c     = sqrt(E/rho);
dtime = 1e-6;
time  = 1e-2;
t     = 0;

nsteps = floor(time/dtime);

istep = 1;

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy
%%
while ( t < time )
    disp(['time step ',num2str(t)])
    
    nvelo(:)     = 0;
    nmassS(:)     = 0;
    nmomentumS(:) = 0;
    % loop over bodies (update nodal momenta without contact)
    for ib=1:1
        %% reset grid data (body contribution)
        nmass(:)     = 0;
        nmomentum(:) = 0;
        niforce(:)   = 0;
        
        body      = bodies(ib);
        elems     = body.elements;
        mpoints   = body.mpoints;
        for ie=1:length(elems)         % loop over computational cells or elements
            e     = elems(ie);
            esctr = element(e,:);      % element connectivity
            enode = node(esctr,:);     % element node coords
            mpts  = mpoints{e};        % particles inside element e
            for p=1:length(mpts)       % loop over particles
                pid    = mpts(p);
                xp     = body.coord(pid,:);
                Mp     = body.mass(pid);
                vp     = body.velo(pid,:);
                Vp   = body.volume(pid);
                stress =  bodies(ib).stress(pid,:);
                
                pt(1)= (2*xp(1)-(enode(1,1)+enode(2,1)))/deltax;
                pt(2)= (2*xp(2)-(enode(2,2)+enode(3,2)))/deltay;
                [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
                J0       = enode'*dNdxi;             % element Jacobian matrix
                invJ0    = inv(J0);
                dNdx     = dNdxi*invJ0;
                % loop over nodes of current element "ie"
                for i=1:length(esctr)
                    id    = esctr(i);
                    dNIdx = dNdx(i,1);
                    dNIdy = dNdx(i,2);
                    nmass(id)       = nmass(id)       + N(i)*Mp;
                    nmomentum(id,:) = nmomentum(id,:) + N(i)*Mp*vp;
                    niforce(id,1)   = niforce(id,1) - Vp*(stress(1)*dNIdx + stress(3)*dNIdy);
                    niforce(id,2)   = niforce(id,2) - Vp*(stress(3)*dNIdx + stress(2)*dNIdy);
                end
            end
        end
        
        activeNodes = bodies(ib).nodes;
               
        % store system momentum and mass
        nmomentumS(activeNodes,:) = nmomentumS(activeNodes,:) + nmomentum(activeNodes,:);
        nmassS    (activeNodes  ) = nmassS    (activeNodes  ) + nmass(activeNodes);
        niforceS  (activeNodes,:) = niforceS  (activeNodes,:) + niforce(activeNodes,:);
    end
           
    nmomentumS = nmomentumS + dtime*niforce;
    nmomentumS(left,:)  = 0;
    nmomentumS(right,:) = 0;

    % update particle velocity and map back to grid velocity
    for ib=1:1
        body      = bodies(ib);
        elems     = body.elements;
        mpoints   = body.mpoints;        
        indices   = 2*ib-1:2*ib;                            
        % loop over computational cells or elements
        for ie=1:length(elems)
            e     = elems(ie);
            esctr = element(e,:);      % element connectivity
            enode = node(esctr,:);     % element node coords
            mpts  = mpoints{e};       % particles inside element e
            % loop over particles
            for p=1:length(mpts)
                pid  = mpts(p);
                xp   = body.coord(pid,:);
                Mp   = body.mass(pid);
                vp   = body.velo(pid,:);
                Vp   = body.volume(pid);
                
                pt(1)= (2*xp(1)-(enode(1,1)+enode(2,1)))/deltax;
                pt(2)= (2*xp(2)-(enode(2,2)+enode(3,2)))/deltay;
                
                [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions

                for i=1:length(esctr)
                    id = esctr(i);                                      
                    vp = vp  + N(i)*(nmomentumS(id,:)-nmomentum(id,:))/nmass(id);                                        
                end
                %@@@ THIS IS WRONG???
                for i=1:length(esctr)
                  id = esctr(i);                  
                  nvelo(id,:)  = nvelo(id,:)  + N(i)*Mp*vp;
                end
                
                bodies(ib).velo(pid,:) = vp;                                
            end
        end
    end
    
    nvelo(left,:)  = 0;
    nvelo(right,:) = 0;
    
    % update particle  position and stresses
    k = 0; u = 0;
    for ib=1:1
        body      = bodies(ib);
        elems     = body.elements;
        mpoints   = body.mpoints;        
        indices   = 2*ib-1:2*ib;                            
        % loop over computational cells or elements
        for ie=1:length(elems)
            e     = elems(ie);
            esctr = element(e,:);      % element connectivity
            enode = node(esctr,:);     % element node coords
            mpts  = mpoints{e};       % particles inside element e
            % loop over particles
            for p=1:length(mpts)
                pid  = mpts(p);
                xp   = body.coord(pid,:);
                Mp   = body.mass(pid);                
                Vp   = body.volume(pid);
                
                pt(1)= (2*xp(1)-(enode(1,1)+enode(2,1)))/deltax;
                pt(2)= (2*xp(2)-(enode(2,2)+enode(3,2)))/deltay;
                
                [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
                J0       = enode'*dNdxi;             % element Jacobian matrix
                invJ0    = inv(J0);
                dNdx     = dNdxi*invJ0;
                
                Lp   = zeros(2,2);
                for i=1:length(esctr)
                    id = esctr(i);       vI=[0 0];
                    if (nmass(id)>0)
                    vI  = nvelo(id,:)/nmass(id);
                    %@@@ ATTENTION: mapped back grid velo only for Lp
                    xp  = xp  + dtime * N(i)*nmomentumS(id,:)/nmass(id);
                    end
                    Lp  = Lp + vI'*dNdx(i,:);
                end
                                
                bodies(ib).coord(pid,:)= xp;
                
                % update stress last
         
                F       = ([1 0;0 1] + Lp*dtime)*reshape(bodies(ib).deform(pid,:),2,2);
                bodies(ib).deform(pid,:) = reshape(F,1,4);
                bodies(ib).volume(pid  ) = det(F)*bodies(ib).volume0(pid);
                dEps    = dtime * 0.5 * (Lp+Lp');
                dsigma  = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
                bodies(ib).stress(pid,:)  = bodies(ib).stress(pid,:) + dsigma';
                bodies(ib).strain(pid,:)  = bodies(ib).strain(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
                
                % compute strain, kinetic energies
                vp   = bodies(ib).velo(pid,:);
                k = k + 0.5*(vp(1)^2+vp(2)^2)*Mp;
                u = u + 0.5*Vp*bodies(ib).stress(pid,:)*bodies(ib).strain(pid,:)';
            end
        end
    end
    
    % update the element particle list
    
    for ib=1:1
        body      = bodies(ib);
        elems     = ones(length(body.volume),1);
        
        for ip=1:length(body.volume)
            x = body.coord(ip,1); y = body.coord(ip,2);
            e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
            elems(ip) = e;
        end
        
        bodies(ib).elements = unique(elems);
        bodies(ib).nodes    = unique(element(bodies(ib).elements,:));
        mpoints = cell(elemCount,1);
        for ie=1:elemCount
            id  = find(elems==ie);
            mpoints{ie}=id;
        end
        
        bodies(ib).mpoints  = mpoints;
    end
    
    % store time,velocty for plotting
    
    ta = [ta;t];
    ka = [ka;k];
    sa = [sa;u];
    
    % VTK output
    
    if (  mod(istep-1,interval) == 0 )
        xp = [bodies(1).coord];
        s  = [bodies(1).stress];
        stress = [s sum(s,2)/3];
        data.stress = stress; 
        vtkFile = sprintf('../../results/mpm/diskWall/%s%d',vtkFileName,istep-1);
        VTKParticles(xp,vtkFile,data);
    end
    
    
    % advance to the next time step
    
    t = t + dtime;
    istep = istep + 1;
end

%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

%%


% pvdFile = fopen(strcat('../results/',vtkFileName,'.pvd'), 'wt');
% 
% fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
% fprintf(pvdFile,'<Collection>\n');
% 
% for i = 1:nsteps
%      if (  mod(i,interval) == 0 )
%     vtuFile = sprintf('%s%d%s',vtkFileName,i,'.vtp');
%     fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
%      end
% end
% 
% fprintf(pvdFile,'</Collection>\n');
% fprintf(pvdFile,'</VTKFile>\n');
% 
% fclose(pvdFile);

%ka = ka*1000;
%sa = sa*1000;

figure 
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
plot(ta(1:end),sa(1:end),'r--','LineWidth',2);
plot(ta(1:end),ka(1:end)+sa(1:end),'g-','LineWidth',2.1);
% xlabel('Time')
% ylabel('Energy (x1E-3)')
% legend('kinetic','strain','total')
% %set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
% %axis([0 3 0 3])


% savefile = 'mpm-explicit-2disks-dt001.mat';
% save(savefile,'ta','ka','sa');

disp([num2str(toc),'   DONE '])

xvi_rigid=[0.03839745962156e0, 0.08169872981078e0, -0.16830127018922e0, -0.21160254037844e0];
yvi_rigid=[-0.26650635094611e0, -0.24150635094611e0, 0.19150635094611e0, 0.16650635094611e0]
