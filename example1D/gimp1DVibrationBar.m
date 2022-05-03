% This file implements the Generalized Material Point Method.
% One dimensional problem with material points.
% The grid is one two-noded linear element.
% This example is taken from "Caveats on the Implementation of the Generalized
% Material Point Method", Buzzi et al, CMES 2008.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% May 2014.

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
L = 25;                % length of the bar
rho = 1;               % density

v0     = 0.1;
n      = 1;            % mode number
c      = sqrt(E/rho);
beta1  = (2*n-1)*0.5*(pi/L);
omega1 = beta1*c;

particleGen = 'Gauss';
particleGen = 'Geom';
cpGIMP      = 1;

%%
%  Computational grid: two-noded elements
m     = [3 4 5 6 7];   % choose number of elements for convenrgence study
elemCount     = 20; % number of elements of the background grid
[mesh]=buildGrid1D(L,elemCount,1); % with ghost cells

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

%% particle size, can be updated to do cpGIMP 

lp = zeros(particles.pCount,1);
lp(:) = deltax/ppc;
lp0   = lp;

%% Time loop

tol = 1e-10;

dtime = 0.1*deltax/c;
time  = (2*pi/omega1)*5;
t     = 0;

ta = [];           % time
va = [];           % velocities
xa = [];           % position
ka = [];           % kinetic energy
sa = [];           % strain energy
er = [];           % errors

nmass      = zeros(nodeCount,1);  % nodal mass vector
nmomentum0 = zeros(nodeCount,1);  % nodal momentum vector
niforce    = zeros(nodeCount,1);  % nodal internal force vector
neforce    = zeros(nodeCount,1);  % nodal external force vector

while ( t < time )
    disp(['time step ',num2str(t)])
    nmass(:)     = 0;
    nmomentum0(:)= 0;
    niforce(:)   = 0;
    % loop over computational cells or elements
    for e=1:elemCount
        esctr = gimpElement{e};        
        mpts  = mpoints{e};
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            lpp = lp(pid,:);
            % loop over nodes contribute to this particle
            for i=1:length(esctr)
                id    = esctr(i);
                x     = xp(pid,:) - nodes(id,:);
                % shape functions and first derivatives
                [phi,dphi]    = getGIMP(x,deltax,lpp);
                % particle mass and momentum to node
                nmass(id)     = nmass(id)     + phi*Mp(pid);
                nmomentum0(id)= nmomentum0(id) + phi*Mp(pid)*vp(pid);
                niforce(id)   = niforce(id)   - Vp(pid)*s(pid)*dphi;
            end
        end
    end
    
    % update nodal momenta
    
    % be careful with left ghost cell
    
    %nmomentum0(1)  = -nmomentum0(3); % Boundary conditions
    nmomentum0(2)  = 0; % Boundary conditions node on symmetry line
    
%     niforce(nodeCount)    = 0;
% 
%     nmomentum(nodeCount)  = 0; % Boundary conditions f1 = m1*a1, a1=0
    
    nmomentum = nmomentum0 + niforce*dtime;
    
    %nmomentum(1)  = -nmomentum(3); % Boundary conditions
    nmomentum(2)  = 0; % Boundary conditions f1 = m1*a1, a1=0
 
    % update particle velocity and position and stresses
    k = 0;
    u = 0;
    for e=1:elemCount
        esctr = gimpElement{e};
        mpts  = mpoints{e};
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            Lp = 0;      
            xxp = xp(pid,:);
            lpp = lp(pid,:);
            for i=1:length(esctr)
                id    = esctr(i);
                x     = xxp - nodes(id,:);
                % shape functions and first derivatives
                [phi,dphi]=getGIMP(x,deltax,lpp);
               
                if nmass(id) > 0
                    vI      = nmomentum(id)/nmass(id);  
                    if ( id == 1 ), vI = - nmomentum(3)/nmass(3); end
                    vp(pid) = vp(pid) +  phi*(nmomentum(id)-nmomentum0(id))/nmass(id);                                        
                    xp(pid) = xp(pid) + dtime * phi*vI;
                else
                    vI = 0.;
                end                                                                
                % gradient velocity                
                Lp = Lp + dphi * vI;
            end
            Fp(pid) = (1 + Lp*dtime)*Fp(pid);
            Vp(pid) = Fp(pid)*Vp0(pid);
            if ( cpGIMP ), lp(pid) = lp0(pid) * Fp(pid); end
            dEps    = dtime * Lp;
            s(pid)  = s(pid)   + E * dEps;
            eps(pid)= eps(pid) + dEps;
            
            k = k + 0.5*vp(pid)^2*Mp(pid);
            u = u + 0.5*s(pid)*eps(pid)*Vp(pid);
        end
    end

    
    % update the element particle list    
    pe = floor(xp/deltax)+1;
    
    for e=1:elemCount
        id  = find(pe==e);
        mpoints{e}=id;
    end

    % advance to the next time step
    t = t + dtime;
    cv = 1/sum(Mp)*(dot(Mp,vp));
    ta = [ta;t];
    va = [va;cv];
    ka = [ka;k];
    sa = [sa;u];
    
    sum(nmass)
    sum(Mp)
end

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
% max(err)
% G=0.0001
disp=[1.193822692497398e-01;
    3.978479481762531e-04;
    4.197781030711051e-04;
    1.478400791855386e-07;
    4.064024736697518e-08];


size=[0.125000000000000;
    0.062500000000000;
    0.031250000000000;
    0.015625000000000;
    0.007812500000000];

polyfit(log(size),log(disp),1)

figure
loglog(size,disp,'black*-','LineWidth',1.8)
hold on
xlabel('Element size')
ylabel('Error')
%legend('MPM','Exact')
set(gca,'FontSize',16)
grid on
legend('G=0.0001','G=0.01')
%set(gca ,'YTickLabel',num2str(disp,1))
% axis([0 100 -0.15 0.2])

%%
figure
plot(ta,er);
set(gca,'FontSize',16)