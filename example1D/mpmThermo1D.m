% This file implements the Material Point Method of Sulsky 1994.
% One dimensional problem with material points.
% The grid is one two-noded linear element.
% Leapfrog time integration.
%
% Thermo MPM 1D example from Tao et al, 2016.
%
% Vinh Phu Nguyen
% Monash University, Australia
% 25 September 2019.

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
clear

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

tic;
%%
%
c      = 1;              % specific heat
rho    = 1;              % density
k      = 1;              % thermal conductivity
L      = 1;              % problem dimension
T1p    = 1;              % prescribed temperature
T0p    = 0;              % initial temperature

% blending PIC/FLIP parameter: 1(FLIP), 0(PIC)
alpha         = 1;
doubleMapping = 1;
beta          = 1; % 1: standard particle position update
% 0: Leroch particle position update



ne          = 100;
[mesh]      = buildGrid1D(L,ne,0);

elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
elements  = mesh.element;
nodes     = mesh.node;
deltax    = mesh.deltax;
deltaxI   = 1 / deltax;
dN1       = -deltaxI;
dN2       =  deltaxI;
%%
% Material points:
ppc       = 1; % # particles per cell
particles = buildParticlesGeom(mesh,ppc,rho);

pCount = particles.pCount;
xp   = particles.xp;        % updated position
xp0  = particles.xp;        % original position
T    = particles.T;         % temperature
C    = particles.C;         % specific heat c
q    = particles.q;         % heat flux -k nalba T
Vp   = particles.Vp;
Vp0  = particles.Vp0;
Mp   = particles.Mp;
C(:) = c;
T(:) = T0p;

%%
% hold on
% plot(nodes,zeros(nodeCount,1)+1/2,'r-s');
% plot(xp,zeros(particles.pCount,1)+1/2,'b*');
% axis([0 L+0.1 0 1.])

%%
% data structure to store the material points for each element
% this data structure is updated for every time step

pElems  = ones(particles.pCount ,1);
mpoints = cell (elemCount ,1);

for p=1:particles.pCount
    x = xp(p);
    e = floor(x/deltax) + 1;
    pElems(p) = e; % particle "p" stays in element "e"
    for e=1:elemCount
        id = find(pElems==e);
        mpoints{e}=id ; % mpoints{e}?> indices of particles in "e"
    end
end

%% Time loop
tol = 0;

% dtime <= deltax * deltax / k
dtime = 1.5*(deltax*deltax/k);
time  = 0.05;
t     = 0.;
istep = 0;


nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,1);  % nodal momentum vector (final)
nmomentum0= zeros(nodeCount,1);  % nodal momentum vector (begin)
ntempMUSL = zeros(nodeCount,1);  % mapped back nodal temperature
niforce   = zeros(nodeCount,1);  % nodal internal force vector
neforce   = zeros(nodeCount,1);  % nodal external force vector

nsteps = floor(time/dtime);
err    = zeros(nsteps,1);
tempe   = zeros(particles.pCount,nsteps); % temperature at particles for all steps
gtempe  = zeros(nodeCount,nsteps);        % temperature at nodes for all steps

while ( t < time )
    disp(['time step ',num2str(t)]);
    nmass(:)     = 0;
    nmomentum0(:)= 0;
    ntempMUSL(:) = 0;
    niforce(:)   = 0;
    neforce(:)   = 0;
    % loop over computational cells or elements
    for e=1:elemCount
        esctr = elements(e,:);
        enode = nodes(esctr);
        mpts  = mpoints{e};
        node1 = esctr(1);
        node2 = esctr(2);
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            xx   = xp(pid);    % position
            mm   = Mp(pid);    % mass
            cp   = C(pid);     % specific heat
            vol  = Vp(pid);    % volume
            tp   = T(pid);     % temperature
            qp   = q(pid);     % heat flux q = - k grad phi T
            % shape functions and first derivatives (not computing der as
            % they are constant)
            N1  = 1 - abs(xx-enode(1))*deltaxI;
            N2  = 1 - abs(xx-enode(2))*deltaxI;
            % particle mass and momentum to node
            nmass(node1)      = nmass(node1)      + N1*mm*cp;
            nmass(node2)      = nmass(node2)      + N2*mm*cp;
            nmomentum0(node1) = nmomentum0(node1) + N1*mm*cp*tp;
            nmomentum0(node2) = nmomentum0(node2) + N2*mm*cp*tp;
            % internal force
            niforce(node1)    = niforce(node1) + vol*qp*dN1;
            niforce(node2)    = niforce(node2) + vol*qp*dN2;
            % external force due to heat sources
        end
    end
    
    % update nodal momenta
    nmomentum = nmomentum0 + niforce*dtime;
    
    % Dirichlet Boundary conditions
    
    ntemp0 = nmomentum0 ./ nmass;
    ntemp  = nmomentum  ./ nmass;
    
    ntemp0(nodeCount) = T1p;
    ntemp (nodeCount) = T1p;
    
    if ( doubleMapping == 1 )
        % update particle velocity, double mapping to get nodal temperature
        for e=1:elemCount
            esctr = elements(e,:);
            node1 = esctr(1);
            node2 = esctr(2);
            enode = nodes(esctr);
            mpts  = mpoints{e};
            % loop over particles
            for p=1:length(mpts)
                pid  = mpts(p);
                xx   = xp(pid);
                mp   = Mp(pid);
                cp   = C(pid);     % specific heat
                N1   = 1 - abs(xx-enode(1))*deltaxI;
                N2   = 1 - abs(xx-enode(2))*deltaxI;
                
                t1   = ntemp(node1);
                t2   = ntemp(node2);
                t10  = ntemp0(node1);
                t20  = ntemp0(node2);
                
                % particle temperature update
           
                T(pid)            = T(pid) + N1*(t1-t10) + N2*(t2-t20);
                
                ntempMUSL(node1)  = ntempMUSL(node1)  + N1*mp*cp*T(pid);
                ntempMUSL(node2)  = ntempMUSL(node2)  + N2*mp*cp*T(pid);
            end
        end
        ntempMUSL = ntempMUSL./nmass;
        ntempMUSL(nodeCount) = T1p;
    else
        ntempMUSL = ntemp;
    end
    
    % update particle temperature and flux (G2P)
    
    for e=1:elemCount
        esctr = elements(e,:);
        node1 = esctr(1);
        node2 = esctr(2);
        enode = nodes(esctr);
        mpts  = mpoints{e};
        T10   = ntemp0(node1);
        T1    = ntemp (node1);
        T20   = ntemp0(node2);
        T2    = ntemp (node2);
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            xx   = xp(pid);
            N1   = 1 - abs(xx-enode(1))*deltaxI;
            N2   = 1 - abs(xx-enode(2))*deltaxI;
            dN1  = -deltaxI;
            dN2  =  deltaxI;
            
            % update particle temperature
            if ( doubleMapping == 0 )
                T(pid)  = T(pid) +  N1*(T1-T10) + N2*(T2-T20);
                q(pid)  = -k * (dN1 * T1 + dN2 * T2);
            else
                q(pid)  = -k * (dN1 * ntempMUSL(node1) + dN2 * ntempMUSL(node2));
            end
        end
    end
    
    % update the element particle list
%     pe = floor(xp/deltax)+1;
%     for e=1:elemCount
%         id  = find(pe==e);
%         mpoints{e}=id;
%     end
    
    % advance to the next time step
    t     = t + dtime;
    istep = istep + 1;
    tempe(:,istep) = T;
    gtempe(:,istep) = ntempMUSL;
end
disp([num2str(toc),'   DONE ']);
%%

% exact solution

exact = zeros(pCount,nsteps);

for p=1:pCount
    xx = xp(p);
    for ti=1:nsteps
        tt = ti*dtime;
        sol = 0;
        for n=1:1000
            alpha =  pi*pi*n*n;
            beta  =  power(-1,n)/(pi*n);
            sol   = sol + beta * exp (-alpha*tt) * sin (n*pi*xx);
        end
        exact(p,ti) = T0p + (T1p-T0p)*xx + 2*(T1p-T0p) * sol;
    end
end

%%
% plot particle temperature at some time instances
% together with the exact solution

set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter',       'latex');
set(0, 'defaultlinelinewidth',2)
set(0, 'DefaultAxesFontSize',20)

figure
set(gca,'FontSize',14)
hold on
plot(xp,tempe(:,50),'black-','LineWidth',2.0);
plot(xp,exact(:,50),'blacks','LineWidth',2.0);
plot(xp,tempe(:,200),'blue-','LineWidth',2.0);
plot(xp,exact(:,200),'blues','LineWidth',2.0);
plot(xp,tempe(:,nsteps),'red-','LineWidth',2.0);
plot(xp,exact(:,nsteps),'reds','LineWidth',2.0);
xlabel('Time')
ylabel('Temperature')
legend('MPM','Exact')
ax=gca;
ax.XAxis.TickLabelFormat='%,.2f';
ax.YAxis.TickLabelFormat='%,.2f';
grid on
% Create smaller axes in top right, and plot on it
% Store handle to axes 2 in ax2.
ax2 = axes('Position',[.24 .4 .4 .4]);
box on;
hold on
range=70:90;
plot(xp(range),tempe(range,50),'black-','LineWidth',2.0);
plot(xp(range),exact(range,50),'blacks','LineWidth',2.0);
grid on;
ax2.XAxis.TickLabelFormat='%,.2f';
ax2.YAxis.TickLabelFormat='%,.2f';

%%
% plot grid temperature at some time instances
% together with the exact solution

set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter',       'latex');
set(0, 'defaultlinelinewidth',2)
set(0, 'DefaultAxesFontSize',20)

figure
set(gca,'FontSize',14)
hold on
plot(nodes,gtempe(:,50),'black-','LineWidth',2.0);
plot(xp,exact(:,50),'blacks','LineWidth',2.0);
plot(nodes,gtempe(:,200),'blue-','LineWidth',2.0);
plot(xp,exact(:,200),'blues','LineWidth',2.0);
plot(nodes,gtempe(:,nsteps),'red-','LineWidth',2.0);
plot(xp,exact(:,nsteps),'reds','LineWidth',2.0);
xlabel('Time')
ylabel('Temperature')
legend('MPM','Exact')
ax=gca;
ax.XAxis.TickLabelFormat='%,.2f';
ax.YAxis.TickLabelFormat='%,.2f';
set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1.0])
grid on
% Create smaller axes in top right, and plot on it
% Store handle to axes 2 in ax2.
ax2 = axes('Position',[.24 .45 .4 .4]);
box on;
hold on
range=50:80;
plot(nodes(range),gtempe(range,50),'black-','LineWidth',2.0);
plot(xp(range),exact(range,50),'blacks','LineWidth',2.0);
grid on;
ax2.XAxis.TickLabelFormat='%,.2f';
ax2.YAxis.TickLabelFormat='%,.2f';