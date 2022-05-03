%

%%

addpath ../grid/
addpath ../basis/
addpath ../particleGen/
addpath ../constitutiveModels/
addpath ../util/
addpath ../geoMesh/
addpath ../externals/
addpath ../postProcessing/


%%
clc
clear all

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%%
%
E = 100;               % Young modulus
L = 6;                % length of the bar
rho = 1;               % density

v0     = 0.1;
n      = 1;            % mode number
c      = sqrt(E/rho);
beta1  = (2*n-1)*0.5*(pi/L);
omega1 = beta1*c;

particleGen = 'Gauss';
particleGen = 'Geom';

%%
%  Computational grid: two-noded elements

elemCount     = 3; % number of elements of the background grid
[mesh]        = buildGrid1D(L,elemCount,1); % with ghost cells

elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
elements  = mesh.element;
nodes     = mesh.node;
deltax    = mesh.deltax;

%%
% Material points:
ppc = 1; % # particles per cell
% switch particleGen
%     case 'Gauss'
%         particles=buildParticles(mesh,ppc,rho);
%     case 'Geom'
%         particles=buildParticlesGeom(mesh,ppc,rho);
% end

particles=buildParticlesGeom(mesh,ppc,rho);
xp  = particles.xp;
vp  = particles.vp;
Vp  = particles.Vp;
Vp0 = particles.Vp0;
Fp  = particles.Fp;
s   = particles.s;
eps = particles.eps;
Mp  = particles.Mp;

% initial velocities

for p=1:particles.pCount
    vp(p) = v0*sin(beta1*xp(p));
end

%%
hold on
plot(nodes,zeros(nodeCount,1)+1/4,'r-s','MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',9,'LineWidth',1.1);
plot(particles.xp,zeros(particles.pCount,1)+1/4,'bo','MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',9,'LineWidth',1.1);
%axis([-10 35 0 1/2])


%%
% data structure to store the material points for each element
% this data structure is updated for every time step

pElems  = ones(particles.pCount ,1);
mpoints = cell (elemCount ,1);

for p=1:particles.pCount
    x = xp(p);
    e = point2ElemIndex1D(x,mesh);
    pElems(p) = e; % particle "p" stays in element "e"
    for e=1:elemCount
        id = find(pElems==e);
        mpoints{e}=id ; % mpoints{e}?> indices of particles in "e"
    end
end

%% find GIMP element connectivity

gimpElement = cell(elemCount,1);

for e=1:elemCount
    neighbors      = getNeighbors1D(e, elemCount);
    neighborNodes  = elements(neighbors,:);
    gimpElement{e} = unique(neighborNodes);
end

%% particle size

lp = deltax/ppc;

%% Time loop

tol = 1e-5;

dtime = 0.1*deltax/c;
time  = (2*pi/omega1)*5;
t     = 0;

ta = [];           % time
va = [];           % velocities
xa = [];           % position
ka = [];           % kinetic energy
sa = [];           % strain energy
er = [];           % errors

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,1);  % nodal momentum vector
niforce   = zeros(nodeCount,1);  % nodal internal force vector
neforce   = zeros(nodeCount,1);  % nodal external force vector


disp(['time step ',num2str(t)])
% store time,velocty for plotting
cv  = 1/sum(Mp)*(dot(Mp,vp));
cvE = v0/(beta1*L)*cos(omega1*t); % center of mass velocity
ta = [ta;t];
va = [va;cv];
%ka = [ka;k];
%sa = [sa;u];
er = [er;abs(cv-cvE)];

nmass(:)     = 0;
nmomentum(:) = 0;
niforce(:)   = 0;
% loop over computational cells or elements
for e=1:elemCount
    esctr = gimpElement{e};
    mpts  = mpoints{e};
    % loop over particles
    for p=1:length(mpts)
        pid  = mpts(p);
        % loop over nodes contribute to this particle
        for i=1:length(esctr)
            id    = esctr(i);
            x     = xp(pid,:) - nodes(id,:);
            % shape functions and first derivatives
            [phi,dphi]    = getGIMP(x,deltax,lp);
            % particle mass and momentum to node
            nmass(id)     = nmass(id)     + phi*Mp(pid);
            nmomentum(id) = nmomentum(id) + phi*Mp(pid)*vp(pid);
        end
    end
end
nodelVelo = nmomentum ./ nmass;
%%

% exact solution

%vExact = v0*cos(omega1.*ta)*sin(beta1*0.5*L);
vExact = v0/(beta1*L)*cos(omega1.*ta); % center of mass velocity
uExact = (v0/omega1)*sin(omega1.*ta)*sin(beta1*0.5*L);

% % rho    = 1;
% w      = 1/L*sqrt(E/rho);
% vExact = 0.1*cos(w.*ta);
% xExact = 0.5*exp((0.1/L/w)*sin(w.*ta));

figure
set(gca,'FontSize',14)
hold on
plot(ta,va,'b-','LineWidth',1.6);
plot(ta,vExact,'r--','LineWidth',2);
xlabel('Time')
ylabel('Velocity')
legend('MPM','Exact')
%axis([0 100 -0.1 0.1])

% figure
% set(gca,'FontSize',14)
% hold on
% plot(ta,xa-0.5*L,'b-','LineWidth',1.6);
% plot(ta,uExact,'r--','LineWidth',2);
% xlabel('Time')
% ylabel('Displacement')
% legend('MPM','Exact')
% set(gca,'FontSize',14)
% axis([0 100 -0.15 0.2])
%%
figure
set(gca,'FontSize',14)
hold on
plot(ta,ka,'b-','LineWidth',1.6);
plot(ta,sa,'r--','LineWidth',2);
plot(ta,ka+sa,'g-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy')
legend('kinetic','strain','total')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
%axis([0 100 0 0.08])

%% convergence plot


%polyfit(log(size),log(disp),1)
% 
% figure
% loglog(size,disp,'black*-','LineWidth',1.8)
% hold on
% xlabel('Element size')
% ylabel('Error')
% %legend('MPM','Exact')
% set(gca,'FontSize',16)
% grid on
% legend('G=0.0001','G=0.01')
%set(gca ,'YTickLabel',num2str(disp,1))
% axis([0 100 -0.15 0.2])

%%
