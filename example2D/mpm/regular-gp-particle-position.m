ghostCell = 0;
lx        = 1;
ly        = 1;
noX0      = 1;      % number of elements along X direction
noY0      = 1;      % number of elements along Y direction
[mesh]    = buildGrid2D(lx,ly,noX0,noY0, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;
omegaC    = mesh.deltax*mesh.deltay;

%% generate material points
ppc           = [2 2];
square.x      = [0 1]; % solid is a square 1 x 1
square.y      = [0 1];
[pmesh]       = buildGrid2D(lx,ly,noX0,noX0, ghostCell); 
[res]         = generateMPForRectangle(square,ppc,pmesh);

%res.position(:,1) = res.position(:,1) + 2*mesh.deltax;
%res.position(:,2) = res.position(:,2) + 2*mesh.deltay;

pCount  = size(res.position,1);
volume  = res.volume;
volume0 = res.volume;
mass    = res.volume*rho;
coords  = res.position;
coords0 = coords;
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stress
velo    = zeros(pCount,2);                % velocity
density = rho*ones(pCount,1);             % density
% thermal part
temp0   = zeros(pCount,1);                % temperature (old)
temp    = zeros(pCount,1);                % temperature
CC      = zeros(pCount,1);                % specifit heat
Q       = zeros(pCount,2);                % heat flux q

% initial temperature = 0, so do nothing
CC(:)   = c;
temp(:) = T0;

% initial velocities, initial stress=0
% for p=1:pCount
%   velo(p,1) = ...;
%   velo(p,2) = ...;
% end


%% plot mesh, particles
figure(1)
hold on
plot_mesh(node,element,'Q4','k-',2.);
%plot_mesh(node,element(bodies{1}.elements,:),'Q4','cy-',2.1);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','r-',2.1);
plot(coords(:,1),coords(:,2),'r.','markersize',40);
pgon1 = polyshape([0 0 0.5 0.5],[0.5 0 0 0.5]);
pgon2 = polyshape([0.5 0.5 1 1],[0.5 0 0 0.5]);
pgon3 = polyshape([0.0 0 0.5 0.5],[1 0.5 0.5 1]);
pgon4 = polyshape([0.5 0.5 1 1],[1 0.5 0.5 1]);
%plot(coords(idx,1),coords(idx,2),'blacks','markersize',20);
%plot(node(idx,1),node(idx,2),'bs','markersize',20);
plot(pgon1)
plot(pgon2)
plot(pgon3)
plot(pgon4)
axis off
print('MPM-regular-particle','-painters','-dpdf','-r1000');

%%
noParticle = 2;

[W,Q]=quadrature( noParticle, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

volume = [];
mass   = [];
coord  = [];


for e=1:mesh.elemCount                 % start of element loop
    sctr = mesh.element(e,:);          %  element scatter vector
    pts  = mesh.node(sctr,:);
    
    for q=1:size(W,1)                           % quadrature loop
        pt=Q(q,:);                              % quadrature point
        wt=W(q);                                % quadrature weight
        [N,dNdxi]=lagrange_basis('Q4',pt);
        J0 = pts'*dNdxi;
        x  = N'*pts;
        
            volume  = [volume;wt*det(J0)];
            mass    = [mass; wt*det(J0)*rho];
            coord   = [coord;x];
        
    end
end

figure
hold on
set(gca,'FontSize',14)
plot_mesh(mesh.node,mesh.element,'Q4','k-',2.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(coord(:,1),coord(:,2),'k.','markersize',40);
plot(coords(:,1),coords(:,2),'r.','markersize',40);
%plot(mesh.node(mesh.lNodes,1),mesh.node(mesh.lNodes,2),'r*','markersize',20);
%plot(mesh.node(mesh.rNodes,1),mesh.node(mesh.rNodes,2),'r*','markersize',20);
%plot(mesh.node(mesh.bNodes,1),mesh.node(mesh.bNodes,2),'r*','markersize',20);
axis off

print('MPM-GP-particle','-painters','-dpdf','-r1000');