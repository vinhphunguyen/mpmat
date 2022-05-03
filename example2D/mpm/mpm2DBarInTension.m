% This file implements the Material Point Method
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% Explicit time integration.
% Use rigid body to impose a moving displacement/velocity boundary condition.
%
% Bar in tension.
% Example taken from Kaul and Sulsky 2004 paper.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% June 2014.

%%

addpath ../fem_util/
addpath ../fem-functions/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
%

E   = 1000;
nu  = 0.0;
rho = 10;
v0  = 0.001;   % imposed velocity
bulk  = E/3/(1-2*nu);
tau   = 1e-3; % c Delta t: relaxation time


stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);

tic;

%% Computational grid

ghostCell=0;
lx     = 0.7;
ly     = 0.5;
numx2  = 7;      % number of elements along X direction
numy2  = 5;      % number of elements along Y direction
[mesh] = buildGrid2D(lx,ly,numx2,numy2, ghostCell);
element= mesh.element;
node   = mesh.node;

%% generate material points

%1. solid domain

l      = 0.3;
numx2  = 6;      % number of elements along X direction
numy2  = 6;      % number of elements along Y direction
[pmesh]= buildGrid2D(l,l,numx2,numy2, 0);

pmesh.node(:,2) = pmesh.node(:,2) + 1*mesh.deltax;

pCount  = pmesh.elemCount;                       % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
coords  = zeros(pCount,2);

for e = 1:pmesh.elemCount
    coord = pmesh.node(pmesh.element(e,:),:);
    a     = 0.5*( coord(1,1)*coord(2,2)  - coord(2,1)*coord(1,2) ...
                + coord(2,1)*coord(3,2)  - coord(3,1)*coord(2,2) ...
                + coord(3,1)*coord(4,2)  - coord(4,1)*coord(3,2) ...
                + coord(4,1)*coord(1,2)  - coord(1,1)*coord(4,2) );
    volume(e)  = a;
    mass(e)    = a*rho;
    coords(e,:) = mean(coord); % center of each element=particle
end

coords0 = coords;

body1.volume  = volume;
body1.volume0 = volume;
body1.mass    = mass;
body1.coord   = coords;
body1.deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress  = zeros(pCount,3);                % stress
body1.strain  = zeros(pCount,3);                % strain
body1.velo    = zeros(pCount,2);                % velocity

%2. rigid body that models the imposed velocity

l     = 0.1;
h     = 0.3;
numx2 = 2;      % number of elements along X direction
numy2 = 6;      % number of elements along Y direction
[rmesh]= buildGrid2D(l,h,numx2,numy2, 0);

rmesh.node(:,1) = rmesh.node(:,1) + 3*mesh.deltax;
rmesh.node(:,2) = rmesh.node(:,2) + 1*mesh.deltay;


pCount  = rmesh.elemCount;                       % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
coords  = zeros(pCount,2);


for e = 1:rmesh.elemCount
    coord = rmesh.node(rmesh.element(e,:),:);
    a     = 0.5*( coord(1,1)*coord(2,2)  - coord(2,1)*coord(1,2) ...
        + coord(2,1)*coord(3,2)  - coord(3,1)*coord(2,2) ...
        + coord(3,1)*coord(4,2)  - coord(4,1)*coord(3,2) ...
        + coord(4,1)*coord(1,2)  - coord(1,1)*coord(4,2) );
    volume(e)  = a;
    mass(e)    = a*rho;
    coords(e,:) = mean(coord); % center of each element=particle
end


pCount = length(volume);

body2.volume = volume;
body2.volume0 = volume;
body2.mass   = mass;
body2.coord  = coords;
body2.deform = repmat([1 0 0 1],pCount,1);     % gradient deformation
body2.stress = zeros(pCount,3);                % stress
body2.strain = zeros(pCount,3);                % strain
body2.velo   = zeros(pCount,2);                % velocity
body2.velo(:,1) = v0;

bodies    = cell(2,1);
bodies{1} = body1;
bodies{2} = body2;
bodyCount = length(bodies);


%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

for ib=1:length(bodies)
    body      = bodies{ib};
    coord     = body.coord;
    elems     = ones(size(coord,1),1);
    
    for ip=1:size(coord,1)
        x = coord(ip,1); y = coord(ip,2);
        e = floor(x/mesh.deltax) + 1 + mesh.numx*floor(y/mesh.deltay);
        elems(ip) = e;
    end
    
    bodies{ib}.elements = unique(elems);
    bodies{ib}.nodes    = unique(element(bodies{ib}.elements,:));
    mpoints = cell(mesh.elemCount,1);
    for ie=1:mesh.elemCount
        id  = find(elems==ie);
        mpoints{ie}=id;
    end
    
    bodies{ib}.mpoints  = mpoints;
end

%% node quantities

nmass     = zeros(mesh.nodeCount,1);  % nodal mass vector
nmomentum = zeros(mesh.nodeCount,2);  % nodal momentum vector
niforce   = zeros(mesh.nodeCount,2);  % nodal internal force vector
neforce   = zeros(mesh.nodeCount,2);  % nodal external force vector
nvelo     = zeros(mesh.nodeCount,2);  % nodal velocities

%% find boundary conditions

eps       = 1e-12;
left      = find(abs(node(:,1))<eps);
fixedBoth = left;

%% plot mesh, particles

figure
set(gca,'FontSize',14)
hold on
plot_mesh(node,element,'Q4','k-',1.7);
xp1 = bodies{1}.coord;
xp2 = bodies{2}.coord;
plot(xp1(:,1),xp1(:,2),'k.','markersize',20);
plot(xp2(:,1),xp2(:,2),'r.','markersize',20);
plot(node(left,1),node(left,2),'b.','markersize',15);
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5])
%axis off


%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

cfl   = 0.1;
c     = sqrt(E/rho);
dtime = cfl*mesh.deltax/c;
time  = 10;
t     = 0;

cc    = sqrt(bulk/rho);
tra   = 0.3/cc;

nsteps = floor(time/dtime);

pos1   = cell(nsteps,1);
pos2   = cell(nsteps,1);

stre   = zeros(nsteps,1); % stress and strain at a monitored particle 
stra   = zeros(nsteps,1);

sigxx   = zeros(nsteps,6); % stress and strain at a monitored particle 
coorx   = zeros(nsteps,6);

istep = 1;

while ( t < time )
    disp(['time step ',num2str(t)])
    %% reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    nvelo(:)     = 0;
    
    %% loop over bodies (update nodal momenta without contact)    
    for ib=1:2
        body      = bodies{ib};
        elems     = body.elements;
        mpoints   = body.mpoints;
        % loop over computational cells or elements
        for ie=1:length(elems)
            e     = elems(ie);
            esctr = element(e,:);      % element connectivity
            enode = node(esctr,:);     % element node coords
            mpts  = mpoints{e};       % particles inside element e            
            % loop over particles
            for p=1:length(mpts)
                pid    = mpts(p);
                xp     = body.coord(pid,:);
                stress = body.stress(pid,:);
                Mp     = body.mass(pid);
                vp     = body.velo(pid,:);
                Vp     = body.volume(pid);
                % loop over nodes of current element "ie"
                for i=1:length(esctr)
                    id    = esctr(i);
                    x     = xp - node(id,:);
                    [N,dNdx]=getMPM2D(x,mesh.deltax,mesh.deltay);
                    
                    dNIdx = dNdx(1);
                    dNIdy = dNdx(2);
                    nmass(id)       = nmass(id)       + N*Mp;
                    nmomentum(id,:) = nmomentum(id,:) + N*Mp*vp;
                    niforce(id,1)   = niforce(id,1) - Vp*(stress(1)*dNIdx + stress(3)*dNIdy);
                    niforce(id,2)   = niforce(id,2) - Vp*(stress(3)*dNIdx + stress(2)*dNIdy);
                end
            end
        end
        
        % update nodal momenta
        nmomentum0 = nmomentum;
        activeNodes=bodies{ib}.nodes;
        nmomentum(activeNodes,:) = nmomentum(activeNodes,:) + niforce(activeNodes,:)*dtime;
        % dynamic relaxation with vicous damping (Sulsky)
        % not necessary for this problem!!!
        nmomentum(activeNodes,1) = nmomentum(activeNodes,1)./(1+tau*nmass(activeNodes));
        nmomentum(activeNodes,2) = nmomentum(activeNodes,2)./(1+tau*nmass(activeNodes));
    end
    
    
    nmomentum(fixedBoth,:) = 0; % fixed boundary conditions
    niforce(fixedBoth,:)   = 0; % fixed boundary conditions       
    niforce(bodies{2}.nodes,:)   = 0; % fixed boundary conditions       
    
    nmomentum(bodies{2}.nodes,1) = nmass(bodies{2}.nodes)*(v0);
    nmomentum(bodies{2}.nodes,2) = 0;
    %acce = nmomentum-nmomentum0;
    %acce(fixedBoth,:)   = 0; % fixed boundary conditions  
    %% update particle velocity 
    
%     for ib=1:1
%         body      = bodies{ib};
%         elems     = body.elements;
%         mpoints   = body.mpoints;
%         % loop over computational cells or elements
%         for ie=1:length(elems)
%             e     = elems(ie);
%             esctr = element(e,:);      % element connectivity
%             enode = node(esctr,:);     % element node coords
%             mpts  = mpoints{e};       % particles inside element e
%             % loop over particles
%             for p=1:length(mpts)
%                 pid  = mpts(p);
%                 xp   = body.coord(pid,:);
%                 Mp   = body.mass(pid);
%                 vp   = body.velo(pid,:);
%                 Vp   = body.volume(pid);              
%                 for i=1:length(esctr)
%                     id = esctr(i);
%                     x     = xp - node(id,:);
%                     [N,dNdx]=getMPM2D(x,mesh.deltax,mesh.deltay);
%                     vp  = vp + dtime * N*niforce(id,:)/nmass(id);
%                     %vp  = vp +  N*acce(id,:)/nmass(id);
%                 end
%                
%                 
%                 bodies{ib}.velo(pid,:) = vp;
%                 
%                 for i=1:length(esctr)
%                     id = esctr(i);
%                     x     = xp - node(id,:);
%                     [N,dNdx]=getMPM2D(x,mesh.deltax,mesh.deltay);
%                     nvelo(id,1:2)  = nvelo(id,1:2) + N*Mp*vp;
%                 end   
%             end
%         end
%     end
%     
%     activeNodes = [bodies{1}.nodes; bodies{2}.nodes];
%     
%     nvelo(activeNodes,1) = nvelo(activeNodes,1) ./ nmass(activeNodes);
%     nvelo(activeNodes,2) = nvelo(activeNodes,2) ./ nmass(activeNodes);
%     
%     nvelo(fixedBoth,:)       = 0;        
%     nvelo(bodies{2}.nodes,1) = v0;
%     nvelo(bodies{2}.nodes,2) = 0;
%        
    for ib=1:1
        body      = bodies{ib};
        elems     = body.elements;
        mpoints   = body.mpoints;
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
                xp0  = body.coord(pid,:);
                vp   = body.velo(pid,:);
                Mp   = body.mass(pid);
                Vp   = body.volume(pid);
                
                Lp   = zeros(2,2);
                for i=1:length(esctr)
                    id = esctr(i);
                    x     = xp0 - node(id,:);
                    [N,dNdx]=getMPM2D(x,mesh.deltax,mesh.deltay);
                    %vI  = nvelo(id,1:2);
                    vI  = nmomentum(id,:)/nmass(id);
                    xp  = xp + dtime * N*nmomentum(id,:)/nmass(id); 
                    vp  = vp + dtime * N*niforce(id,:)/nmass(id);
                    Lp  = Lp  + vI'*dNdx;         % particle gradient velocity
                end
                
                
                bodies{ib}.coord(pid,:)= xp;
                bodies{ib}.velo(pid,:) = vp;
                
                F       = ([1 0;0 1] + Lp*dtime)*reshape(bodies{ib}.deform(pid,:),2,2);
                bodies{ib}.deform(pid,:) = reshape(F,1,4);
                bodies{ib}.volume(pid  ) = det(F)*bodies{ib}.volume0(pid);
                dEps    = dtime * 0.5 * (Lp+Lp');
                bodies{ib}.strain(pid,:)  = bodies{ib}.strain(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
                
                if ( ib == 2 ) % elastic
                    %                     dsigma  = C2 * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
                    %                     bodies{ib}.stress(pid,:)  = bodies{ib}.stress(pid,:) + dsigma';
                else           % plastic
                    dsigma  = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
                    sigma  = bodies{ib}.stress(pid,:) + dsigma';
                    bodies{ib}.stress (pid,:) = sigma;
                end                                                         
            end
        end
    end
    
    bodies{2}.coord = bodies{2}.coord + dtime* bodies{2}.velo;
    
    % update the element particle list
    
    for ib=1:length(bodies)
        body      = bodies{ib};
        coord     = body.coord;
        elems     = ones(size(coord,1),1);
        
        for ip=1:length(elems)
            x = coord(ip,1); y = coord(ip,2);
            e = floor(x/mesh.deltax) + 1 + mesh.numx*floor(y/mesh.deltay);
            elems(ip) = e;
        end
        
        bodies{ib}.elements = unique(elems);
        bodies{ib}.nodes    = unique(element(bodies{ib}.elements,:));
        
        mpoints = cell(mesh.elemCount,1);
        for ie=1:mesh.elemCount
            id  = find(elems==ie);
            mpoints{ie}=id;
        end
        
        bodies{ib}.mpoints  = mpoints;
    end
    
    pos1{istep} = bodies{1}.coord;
    pos2{istep} = bodies{2}.coord;
    
    % record uniaxial stress and strain
    stre(istep) = bodies{1}.stress(6,1);
    stra(istep) = bodies{1}.strain(6,1);
    
    sigxx(istep,:) = bodies{1}.stress([1:6],1)';
    coorx(istep,:) = bodies{1}.coord ([1:6],1);
    
    % advance to the next time step
    
    t = t + dtime;
    istep = istep + 1;
    
    nvelo0 = nvelo;
end

disp([num2str(toc),'   POST-PROCESSING '])

%%

step = nsteps;
figure
set(gca,'FontSize',14)
hold on
plot_mesh(node,element,'Q4','k-',1.7);
xp1 = pos1{step};
xp2 = pos2{step};
plot(xp1(:,1),xp1(:,2),'k.','markersize',20);
plot(xp2(:,1),xp2(:,2),'r.','markersize',20);
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5])
%axis off

%%

figure
set(gca,'FontSize',14)
hold on
plot(stra*1000,stre,'b-','LineWidth',1.6);
%plot(ss.ta(1:end),ss.ka(1:end),'black-','LineWidth',2.1);
xlabel('strain (x 0.001)')
ylabel('stress')
%legend('kinetic','strain','total')
set(gca,'XTick',[0 5 10 15 20 25 30 35])
set(gca,'YTick',[0 5 10 15 20 25 30 35])
%axis([0 3 0 3])
axis equal

%%

figure
set(gca,'FontSize',14)
hold on
plot(coorx(1000,:),sigxx(1000,:),'black-','LineWidth',1.6);
plot(coorx(2000,:),sigxx(2000,:),'black--','LineWidth',1.6);
plot(coorx(3000,:),sigxx(3000,:),'red-','LineWidth',1.6);
plot(coorx(4000,:),sigxx(4000,:),'red--','LineWidth',1.6);
plot(coorx(5000,:),sigxx(5000,:),'b-','LineWidth',1.6);
plot(coorx(6000,:),sigxx(6000,:),'b--','LineWidth',1.6);
xlabel('x')
ylabel('stress')
legend('1s','2s','3s','4s','5s','6s')
set(gca,'XTick',[0 5 10 15 20 25 30 35])
set(gca,'YTick',[0 5 10 15 20 25 30 35])
%axis([0 3 0 3])
axis equal

disp([num2str(toc),'   DONE '])
