% This file implements the Material Point Method of Sulsky 1994.
% One dimensional problem with material points.
% The grid is one two-noded linear element.
% This example is taken from "Caveats on the Implementation of the Generalized
% Material Point Method", Buzzi et al, CMES 2008.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% June 2013.
% Modified:
% Phu Nguyen, Monash Uni, 11 September 2019
% Add: blending PIC/FLIP for particle velocity/position update
% Add: matlab2tikz for better integration with LaTEX

%%

addpath ../grid/
addpath ../basis/
addpath ../particleGen/
addpath ../constitutiveModels/
addpath ../util/
addpath ../geoMesh/
addpath ../externals/
addpath ../postProcessing/
addpath('/Users/vingu/my-codes/matlab2tikz/src');
%matlab2tikz('sampleFigure.tex','height','0.5\textwidth');
opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)
%%
clc
clear all

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%%
%
E   = 100;               % Young modulus
L   = 25;                % length of the bar
rho = 1;                 % density
v0  = 0.1;
n   = 1;                 % mode number

c      = sqrt(E/rho);
beta1  = (2*n-1)*0.5*(pi/L);
omega1 = beta1*c;

% blending PIC/FLIP parameter: 1(FLIP), 0(PIC)
alpha         = 1;
doubleMapping = 0;
beta          = 1; % 1: standard particle position update
                   % 0: Leroch particle position update
%%
%  Computational grid: two-noded elements

noElems   = 16;
[mesh]    = buildGrid1D(L,noElems,0);

elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
elements  = mesh.element;
nodes     = mesh.node;
deltax    = mesh.deltax;
deltaxI   = 1/deltax;

%%
% Material points:
ppc       = 2; % # particles per cell
particles = buildParticlesGeom(mesh,ppc,rho);

xp  = particles.xp;      % position
vp  = particles.vp;      % velocity
Vp  = particles.Vp;      % volume
Vp0 = particles.Vp0;     % intial volume
Fp  = particles.Fp;      % deformation gradient
s   = particles.s;       % stress
eps = particles.eps;     % strain
Mp  = particles.Mp;      % mass

% initial velocities

for p=1:particles.pCount
    vp(p) = v0*sin(beta1*xp(p));
end
vp0=vp;
%q  = Mp*vp;                          % momentum

%%
hold on
plot(nodes,zeros(nodeCount,1)+1/2,'r-s');
plot(particles.xp,zeros(particles.pCount,1)+1/2,'b*');
axis([0 L+0.5 0 1.])

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

tol = 1e-5;

dtime = 0.2*deltax/c;
time  = (2*pi/omega1)*5; %time=10*dtime;
t     = 0;

ta = [];           % time
va = [];           % velocities
xa = [];           % position
ka = [];           % kinetic energy
sa = [];           % strain energy

nmass      = zeros(nodeCount,1);  % nodal mass vector at time 't'
nmomentum0 = zeros(nodeCount,1);  % nodal momentum vector at time 't'
nmomentum  = zeros(nodeCount,1);  % nodal momentum vector at time 't+dt'
niforce    = zeros(nodeCount,1);  % nodal internal force vector at 't'
neforce    = zeros(nodeCount,1);  % nodal external force vector at 't'
nvelo      = zeros(nodeCount,1);  % store double mapped velocity

%% time loop
while ( t < time )
    % reset grid data
    nmass(:)      = 0;
    nmomentum0(:) = 0;
    nvelo(:)      = 0;
    niforce(:)    = 0;
    % Particles to grid (P2G)
    % loop over computational cells or elements
    for e=1:elemCount
        esctr = elements(e,:);
        node1 = esctr(1);
        node2 = esctr(2);
        enode = nodes(esctr);
        mpts  = mpoints{e};
        % loop over particles
        for p=1:length(mpts)
            pid   = mpts(p);
            xx    = xp(pid);
            mp    = Mp(pid);
            vol   = Vp(pid);
            velo  = vp(pid);
            sigma = s(pid);
            % shape functions and first derivatives
            N1  = 1 - abs(xx-enode(1))*deltaxI;
            N2  = 1 - abs(xx-enode(2))*deltaxI;
            dN1 = -deltaxI;
            dN2 =  deltaxI;
            % particle mass and momentum to node
            nmass(node1)      = nmass(node1)     + N1*mp;
            nmass(node2)      = nmass(node2)     + N2*mp;
            nmomentum0(node1) = nmomentum0(node1) + N1*mp*velo;
            nmomentum0(node2) = nmomentum0(node2) + N2*mp*velo;
            % internal force
            niforce(node1) = niforce(node1) - vol*sigma*dN1;
            niforce(node2) = niforce(node2) - vol*sigma*dN2;
        end
    end
    
    % update nodal momenta
    nmomentum = nmomentum0 + niforce*dtime;
    % Dirichlet Boundary conditions
    nmomentum0(1)  = 0;
    nmomentum (1)  = 0;
    
    if ( doubleMapping == 1 )
        % update particle velocity, double mapping to get nodal velocities
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
                N1   = 1 - abs(xx-enode(1))*deltaxI;
                N2   = 1 - abs(xx-enode(2))*deltaxI;
                
                v1   = nmomentum(node1)/nmass(node1);
                v2   = nmomentum(node2)/nmass(node2);
                
                % particle velocity blending PIC/FLIP
                % beta = 0: PIC
                % beta = 1: FLIP (most MPM uses this)
                vp(pid)  = alpha  * ( vp(pid) + N1*(nmomentum(node1)-nmomentum0(node1))/nmass(node1) + ...
                                              + N2*(nmomentum(node2)-nmomentum0(node2))/nmass(node2) ) + ...
                    (1-alpha) * ( N1*v1 + N2*v2 );
                
                nvelo(node1)  = nvelo(node1)  + N1*mp*vp(pid);
                nvelo(node2)  = nvelo(node2)  + N2*mp*vp(pid);
            end
        end
        nvelo = nvelo./nmass;
        nvelo(1) = 0;
    else
        nvelo = nmomentum./nmass;
    end
    
 
    % update particle velocity and position and stresses (G2P)
    k = 0;
    u = 0;
    for e=1:elemCount
        esctr = elements(e,:);
        node1 = esctr(1);
        node2 = esctr(2);
        enode = nodes(esctr);
        mpts  = mpoints{e};
        v10   = nmomentum0(node1)/nmass(node1);
        v1    = nmomentum (node1)/nmass(node1);
        v20   = nmomentum0(node2)/nmass(node2);
        v2    = nmomentum (node2)/nmass(node2);
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            xx   = xp(pid);
            N1   = 1 - abs(xx-enode(1))*deltaxI;
            N2   = 1 - abs(xx-enode(2))*deltaxI;
            dN1  = -deltaxI;
            dN2  =  deltaxI;
            
            % update particle velocity/position
            %vp(pid)  = alpha*vp(pid) +  N1*(v1-alpha*v10) + N2*(v2-alpha*v20);
            
            if ( doubleMapping == 0 ) 
                vp(pid)  = alpha*vp(pid) +  N1*(v1-alpha*v10) + N2*(v2-alpha*v20);
                % gradient velocity
                Lp = dN1 * v1 + dN2 * v2;
            else
                Lp = dN1 * nvelo(node1) + dN2 * nvelo(node2);
            end
            
            xp(pid)  = xp(pid) + (1-beta)*dtime*vp(pid) + beta*dtime*(N1*v1+N2*v2);
            
           
            
            Fp(pid) = (1 + Lp*dtime)*Fp(pid);
            Vp(pid) = Fp(pid)*Vp0(pid);
            dEps    = dtime * Lp;
            s(pid)  = s(pid)   + E * dEps;
            eps(pid)= eps(pid) + dEps;
            
            k = k + 0.5*vp(pid)^2*Mp(pid);
            u = u + 0.5*s(pid)*eps(pid)*Vp(pid);
        end
    end
    
    % store time,velocty for plotting
    
    
    
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
end
%%


%%

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
%set(gca,'FontSize',14)
hold on
plot(ta,va,'b-','LineWidth',1.4);
plot(ta,vExact,'r--','LineWidth',1.4);
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
%set(gca,'FontSize',14)
hold on
plot(ta,ka,'b-','LineWidth',1.2);
plot(ta,sa,'r--','LineWidth',1.2);
plot(ta,ka+sa,'g-','LineWidth',1.2);
xlabel('Time')
ylabel('Energy')
legend('kinetic','strain','total')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
%axis([0 100 0 0.08])

