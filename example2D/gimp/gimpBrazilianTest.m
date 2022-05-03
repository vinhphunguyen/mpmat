%
% MPM modelling of the Brazilian test.
% Paramters taken from the article "Numerical modelling of a Brazilian Disc
% test of layered rocks using the combined finite-discrete element method"
% by Mahabadi et al, 2009, ROCKENG09.
%
%
% Vinh Phu Nguyen
% University of Adelaide, Adelaide, Australia.
% November 2014.

%%

addpath ../fem_util/
addpath ../fem-functions/
addpath ../post-processing/
addpath mex/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

global mat

%% Material properties
%

mat.E   = 50e3;              % Young's modulus
mat.nu  = 0.25;              % Poisson ratio
rho     = 2700e-12;          % density
mat.stressState ='PLANE_STRESS'; % set  to either 'PLANE_STRAIN' or "PLANE_STRESS'
mat.ft = 8.12;
mat.Gf = 0.05;
mat.ao = elasticityMatrix(mat.E,mat.nu,mat.stressState);
mat.penalty = 1e8; % penalty stiffness for cohesive law in compression
mat.ks = 10;

v     = 0;               % initial particle velocity
v0    = 1000;            % imposed velocity of  the load platen

% OPTIONS FOR THE CONSTITUTIVE MODEL

option.implicit = 0;
option.tolerance= 1e-4;
option.tangent  = 1; %
option.iterMax  = 20;

vtkFileName  = 'brazilianTest';


tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell= 0;
lx       = 6*50;
ly       = 6*50;
numx2    = 50;       % number of elements along X direction
numy2    = 50;       % number of elements along Y direction
[mesh]   = buildGrid2D(lx,ly,numx2,numy2, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;
mesh.H    = sqrt ( mesh.deltax * mesh.deltay );

%% particle generation
% body1

ppc    = 3; % # of particle per cell is ppc x ppc
[W,Q]=quadrature(  ppc, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

dx = mesh.deltax/(ppc);
dy = mesh.deltay/(ppc);


    mesh.lpx = dx;
    mesh.lpy = dy;


center = 0.5*[lx ly];
radius = 1.0*50;

volume = []; mass   = []; coord  = [];

for e=1:elemCount                 % start of element loop
    sctr = element(e,:);          %  element scatter vector
    pts  = node(sctr,:);
    x1 = pts(1,:); % first corner of the cell
    for i=1:ppc
        for j=1:ppc
            x(1) = x1(1) + dx*0.5 + (j-1)*dx;
            x(2) = x1(2) + dy*0.5 + (i-1)*dy;
            r  = norm(x-center);
            if ( r-radius < 0 )
                volume  = [volume;dx*dy];
                mass    = [mass; dx*dy*rho];
                coord   = [coord;x];
            end
        end
    end
end

% coord(:,1) = coord(:,1) + mesh.deltax;
% coord(:,2) = coord(:,2) + mesh.deltay;

pCount = length(volume);

% stored in body1 structure

body1.coord   = coord;
body1.volume  = volume;
body1.volume0 = volume;
body1.mass    = mass;
body1.deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress  = zeros(pCount,3);                % stress
body1.strain  = zeros(pCount,3);                % strain
body1.velo    = ones(pCount,2)*v;               % velocity
body1.C       = mat.ao;
body1.color   = ones(pCount,1);
%body1.deform0 = body1.deform;
body1.gravity = 0;

body1.status   = zeros(1,pCount); % cracked or not at all Gauss points
body1.loading  = zeros(1,pCount); % loading status at all Gauss points
body1.ai       = zeros(4,pCount); % tangent stiffness of localisation band
body1.normal   = zeros(2,pCount); % normal vectors at all Gauss points
body1.kappa    = zeros(1,pCount); % history variables (updated)
body1.jump     = zeros(2,pCount); % jump
body1.traction = zeros(2,pCount); % traction vector (updated), for visualisation
body1.damage   = zeros(1,pCount); % damage variable for visualisation
body1.cracks   = zeros(4,pCount);


%% body 2 and body 3 are the rigid plate for moving boundary conditions
% body 2

volume = []; mass   = []; coord  = [];

for e=1:elemCount                 % start of element loop
    sctr = element(e,:);            %  element scatter vector
    pts  = node(sctr,:);
    
    x1 = pts(1,:); % first corner of the cell
    for i=1:ppc
        for j=1:ppc
            x(1) = x1(1) + dx*0.5 + (j-1)*dx;
            x(2) = x1(2) + dy*0.5 + (i-1)*dy;
            if ( ( x(1) > 0.5*lx -0.5*mesh.deltax ) && ( x(1) < 0.5*lx+0.5*mesh.deltax  ) ...
                    && ( x(2) < 0.5*ly - 50 ) && ( x(2) > 0.5*ly - 50 - mesh.deltay  ) )
                volume  = [volume;dx*dy];
                mass    = [mass; dx*dy*rho];
                coord   = [coord;x];
            end
        end
    end
end

% coord(:,1) = coord(:,1) + mesh.deltax;
% coord(:,2) = coord(:,2) + mesh.deltay;

pCount = length(volume);

body2.volume  = volume;
body2.volume0 = volume;
body2.mass    = mass;
body2.coord   = coord;
body2.deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
body2.stress  = zeros(pCount,3);                % stress
body2.strain  = zeros(pCount,3);                % strain
body2.velo    = ones(pCount,2)*v0;              % velocity
body2.velo(:,1) = 0;
%body2.deform0 = body2.deform;
body2.gravity = 0;
body2.color   = 2*ones(pCount,1);

% body 3

volume = []; mass   = []; coord  = [];

for e=1:elemCount                 % start of element loop
    sctr = element(e,:);            %  element scatter vector
    pts  = node(sctr,:);
    
    x1 = pts(1,:); % first corner of the cell
    for i=1:ppc
        for j=1:ppc
            x(1) = x1(1) + dx*0.5 + (j-1)*dx;
            x(2) = x1(2) + dy*0.5 + (i-1)*dy;
            if ( ( x(1) > 0.5*lx -0.5*mesh.deltax ) && ( x(1) < 0.5*lx+0.5*mesh.deltax  ) ...
                    && ( x(2) > 0.5*ly + 50 ) && ( x(2) < 0.5*ly + 50 + mesh.deltay  ) )
                volume  = [volume;dx*dy];
                mass    = [mass; dx*dy*rho];
                coord   = [coord;x];
            end
        end
    end
end

% coord(:,1) = coord(:,1) + mesh.deltax;
% coord(:,2) = coord(:,2) + mesh.deltay;

pCount = length(volume);

body3.volume  = volume;
body3.volume0 = volume;
body3.mass    = mass;
body3.coord   = coord;
body3.deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
body3.stress  = zeros(pCount,3);                % stress
body3.strain  = zeros(pCount,3);                % strain
body3.velo    = -ones(pCount,2)*v0;             % velocity
body3.velo(:,1) = 0;
%body2.deform0 = body2.deform;
body3.gravity = 0;
body3.color   = 2*ones(pCount,1);

bodies    = cell(3,1);
bodies{1} = body1;
bodies{2} = body2;
bodies{3} = body3;
bodyCount = length(bodies);

%% find GIMP element connectivity

gimpElement = cell(elemCount,1);

for e=1:elemCount
    neighbors      = getNeighbors(e, mesh.numx, mesh.numy);
    neighborNodes  = element(neighbors,:);
    gimpElement{e} = unique(neighborNodes);
end

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

% The folowing function is only for MPM
%bodies = findActiveElemsAndNodes(bodies,mesh);

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
    bodies{ib}.nodes    = unique(vertcat(gimpElement{bodies{ib}.elements}));
    
    mpoints = cell(elemCount,1);
    for ie=1:elemCount
        id  = find(elems==ie);
        mpoints{ie}=id;
    end
    
    bodies{ib}.mpoints  = mpoints;
end


%% node quantities

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,2);  % nodal momentum vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector (no need for this exam)
nvelo     = zeros(nodeCount,2);  % nodal velocity vector
nacce     = zeros(nodeCount,2);  % nodal acceleration vector

%% plot mesh, particles

figure(1)
set(gca,'FontSize',14)
hold on
plot_mesh(node,element,'Q4','k-',1.6);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','cy-',1.6);
xp1 = bodies{1}.coord;
xp2 = bodies{2}.coord;
xp3 = bodies{3}.coord;
plot(xp1(:,1),xp1(:,2),'k.','markersize',15);
plot(xp2(:,1),xp2(:,2),'r.','markersize',15);
plot(xp3(:,1),xp3(:,2),'c.','markersize',15);
plot(node(bodies{1}.nodes,1),node(bodies{1}.nodes,2),'*');
plot(node(bodies{2}.nodes,1),node(bodies{2}.nodes,2),'*');

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance
courant=0.02;

c     = sqrt(mat.E/rho);
dtime = courant*(mesh.deltax/c);
time  = 3.5e-3;% explicit, high velocity
t     = 0;

nsteps   = floor(time/dtime);
interval = floor (nsteps / 50);  % interval for VTK output


save_fre = zeros(round(nsteps/interval),1);
save_dis = zeros(round(nsteps/interval),1);
sstep    = 1;

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

figure(2)

while ( t < time )
    disp(['time step ',num2str(t)])
    %% reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    %nvelo(:)     = 0;
    for ib=1:1                         %% loop over bodies
        body      = bodies{ib};
        elems     = body.elements;
        mpoints   = body.mpoints;
        for ie=1:length(elems)         % loop over computational cells or elements
            e     = elems(ie);
            esctr = gimpElement{e};    % element connectivity
            mpts  = mpoints{e};        % particles inside element e
            for p=1:length(mpts)       % loop over particles
                pid    = mpts(p);
                xp     = body.coord(pid,:);
                stress = body.stress(pid,:);
                Mp     = body.mass(pid);
                vp     = body.velo(pid,:);
                Vp     = body.volume(pid);
                for i=1:length(esctr)  % loop over nodes of current element "ie"
                    id    = esctr(i);
                    x     = xp - node(id,:);
                    [N,dNdx]=getGIMP2D(x,mesh.deltax,mesh.deltay,mesh.lpx,mesh.lpy);
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
        
        activeNodes=bodies{ib}.nodes;
        nmomentum(activeNodes,:) = nmomentum(activeNodes,:) + niforce(activeNodes,:)*dtime;
    end
    
    % boundary conditions from rigid particles
    %nmomentum(bodies{2}.nodes,1) = 0;
    nmomentum(bodies{2}.nodes,2) = nmass(bodies{2}.nodes)*( v0);
    nmomentum(bodies{2}.nodes,1) = 0;
    
    nmomentum(bodies{3}.nodes,2) = nmass(bodies{3}.nodes)*( -v0);
    nmomentum(bodies{3}.nodes,1) = 0;
    
    nforce = niforce;
    
    niforce(bodies{2}.nodes,2)   = 0;
    niforce(bodies{3}.nodes,2)   = 0;
    
    for ib=1:1
        body      = bodies{ib};
        elems     = body.elements;
        mpoints   = body.mpoints;
        for ie=1:length(elems)         % loop over computational cells or elements
            e     = elems(ie);
            esctr = gimpElement{e};    % element connectivity
            mpts  = mpoints{e};        % particles inside element e
            for p=1:length(mpts)       % loop over particles
                pid  = mpts(p);
                xp   = body.coord(pid,:);
                xp0  = body.coord(pid,:);
                vp   = body.velo(pid,:);
                Lp   = zeros(2,2);
                for i=1:length(esctr)
                    id  = esctr(i);
                    x   = xp0 - node(id,:);
                    [N,dNdx]=getGIMP2D(x,mesh.deltax,mesh.deltay,mesh.lpx,mesh.lpy);
                    if (nmass(id)~=0)
                        xp  = xp + dtime * N*nmomentum(id,:)/nmass(id);
                        vp  = vp + dtime * N*niforce(id,:)/nmass(id);
                        Lp  = Lp + (nmomentum(id,:)/nmass(id))'*dNdx;% particle gradient velocity
                    end
                end
                
                bodies{ib}.coord(pid,:)= xp;
                bodies{ib}.velo(pid,:) = vp;
                
                F       = ([1 0;0 1] + Lp*dtime)*reshape(bodies{ib}.deform(pid,:),2,2);
                bodies{ib}.deform(pid,:) = reshape(F,1,4);
                bodies{ib}.volume(pid  ) = det(F)*bodies{ib}.volume0(pid);
                dEps    = dtime * 0.5 * (Lp+Lp');
                
                history.H       = mesh.H;
                history.cracked = body.status (  pid);
                history.loading = body.loading(  pid);
                history.normal  = body.normal (:,pid);
                history.kappa0  = body.kappa  (  pid);
                history.jump0   = body.jump   (:,pid);
                history.ai      = body.ai     (:,pid);
                history.sigma0  = body.stress (  pid,:)';
                history.a       = body.ai     (:,pid);
                
                load.eps0 = bodies{ib}.strain(pid,:)';
                load.dEps = [dEps(1,1); dEps(2,2); 2*dEps(1,2)];
                
                out             = updateTwoScaleCohesiveRK(history,load,option);
                
                bodies{ib}.stress(pid,:)  = out.sigma;
                bodies{ib}.strain(pid,:)  = bodies{ib}.strain(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
                
                bodies{ib}.ai      (:,pid)  = reshape(out.K,1,size(out.K,1)^2);
                bodies{ib}.kappa   (  pid)  = out.kappa;
                bodies{ib}.jump    (:,pid)  = out.u;
                bodies{ib}.loading (  pid)  = out.loading;
                bodies{ib}.damage  (  pid)  = out.damage;
                bodies{ib}.traction(:,pid)  = out.trac;
                
            end
        end
    end
    
    % for rigid bodies (load platens)
    bodies{2}.coord = bodies{2}.coord + dtime* bodies{2}.velo;
    bodies{3}.coord = bodies{3}.coord + dtime* bodies{3}.velo;
    
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
        bodies{ib}.nodes    = unique(vertcat(gimpElement{bodies{ib}.elements}));
        %         if ib==2
        %           bodies{ib}.nodes    = unique(element(bodies{ib}.elements,:));
        %         end
        mpoints = cell(elemCount,1);
        for ie=1:elemCount
            id  = find(elems==ie);
            mpoints{ie}=id;
        end
        
        bodies{ib}.mpoints  = mpoints;
    end
    
    % VTK output
    if (  mod(istep-1,interval) == 0 )
        xp = [bodies{1}.coord;bodies{2}.coord;bodies{3}.coord];
        s1 = [bodies{1}.stress zeros(size(bodies{1}.stress,1),1)];
        s2 = [bodies{2}.stress zeros(size(bodies{2}.stress,1),1)];
        s3 = [bodies{3}.stress zeros(size(bodies{3}.stress,1),1)];
        
        %         for i=1:size(bodies{1}.stress,1)
        %             devSig = I_dev*s1(i,1:3)';
        %             s1(i,4)= sqrt(  devSig(1)^2 + devSig(2)^2 + 2*devSig(3)^2  );
        %         end
        data.stress  = [s1;s2;s3];
        data.velo    = [bodies{1}.velo;  bodies{2}.velo;bodies{3}.velo];
        data.color   = [bodies{1}.color; bodies{2}.color;bodies{3}.color];
        data.damage  = [bodies{1}.damage'; zeros(length(body2.mass),1);zeros(length(body2.mass),1)];
        vtkFile = sprintf('../results/gimp/brazilian/%s%d',vtkFileName,istep-1);
        VTKParticles(xp,vtkFile,data);
        pos{istep} = xp;
        
        % calculate and plot load-disp curve
        sstep           = sstep + 1;
        boundary_layer  = intersect(bodies{1}.nodes,bodies{2}.nodes);
        size_bl         = size(boundary_layer,1)/3;
        boundary_layer  = boundary_layer(size_bl+2:size_bl*2-1);
        save_fre(sstep) = sum(nforce(boundary_layer,2));
        save_dis(sstep) = save_dis(sstep-1) + v0 * dtime;
        
        plot(save_dis(1:sstep),save_fre(1:sstep),'r-.','LineWidth',2);
        drawnow;
        
    end
    
    % check failure
    body      = bodies{1};
    for p=1:length(bodies{1}.mass)                          % quadrature loop
        if body.status(p)==1 continue; end
        sigma  = body.stress(p,:);
        sigmai = getPrincipalStress(sigma,mat.nu,mat.stressState);
        yield  = sigmai(1) - mat.ft;
        ft     = sigmai(1);
        if (yield>=0)
            isFailure     = 1;
            bodies{1}.status(p)   = 1;
            bodies{1}.ai(:,  p)   = [-ft^2/mat.Gf 0 0 mat.ks];
            n             = getPrincipalDirection( sigma );
            %n=[1 0];
            bodies{1}.normal(:,p) = n;
            % compute a crack segment, for visualisation purpose
            xp   = body.coord(p,:);
            tem1 = 0.25*mesh.H*n(2);
            tem2 = 0.25*mesh.H*n(1);
            bodies{1}.cracks(:,p) = [xp(1)-tem1 xp(2)-tem2 xp(1)+tem1 xp(2)+tem2];
        end
    end
    
    % advance to the next time step
    t = t + dtime;
    istep = istep + 1;
end

%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

figure
set(gca,'FontSize',14)
hold on
plot_mesh(node,element,'Q4','k-',1.6);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','cy-',1.6);
xp1 = bodies{1}.coord;
xp2 = bodies{2}.coord;
xp3 = bodies{3}.coord;
plot(xp1(:,1),xp1(:,2),'b.','markersize',15);
plot(xp2(:,1),xp2(:,2),'r.','markersize',15);
plot(xp3(:,1),xp3(:,2),'r.','markersize',15);

for p=1:length(body1.mass) 
    loading = bodies{1}.loading(p);
    if ( loading )
       plot(bodies{1}.cracks([1 3],p),bodies{1}.cracks([2 4],p),...
            'red-','LineWidth',2.8)
    else
        plot(bodies{1}.cracks([1 3],p),bodies{1}.cracks([2 4],p),...
             'black-','LineWidth',1.6)
    end
end

disp([num2str(toc),'   DONE '])
