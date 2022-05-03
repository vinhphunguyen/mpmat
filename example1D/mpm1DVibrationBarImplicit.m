% This file implements the Material Point Method of Sulsky 1994.
% One dimensional problem with material points.
% The grid is one two-noded linear element.
% This example is taken from "Caveats on the Implementation of the Generalized
% Material Point Method", Buzzi et al, CMES 2008.
%
% Newmark time integration scheme.
%
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia.
% 10 August 2015.

%%

addpath ../fem_util/


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

% Newmark parameters
beta  = 0.25;
gamma = 0.5;

%%
%  Computational grid: two-noded elements

[mesh]=buildGrid1D(L,14,0);

elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
elements  = mesh.element;
nodes     = mesh.node;
deltax    = mesh.deltax;

%%
% Material points:
ppc = 4; % # particles per cell
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

dis     = zeros(particles.pCount,1);                % displacement
dis0    = zeros(particles.pCount,1);                % old displacement
vp0     = vp;
acce    = zeros(particles.pCount,1);                % acceleration
acce0   = zeros(particles.pCount,1);                % old acceleration

%q  = Mp*vp;                          % momentum

%%
% hold on
% plot(nodes,zeros(nodeCount,1)+1/2,'r-s');
% plot(particles.xp,zeros(particles.pCount,1)+1/2,'b*');
% axis([0 L+0.5 0 1.])

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

dtime = 0.1*deltax/c;
time  = (2*pi/omega1)*5; time =10*dtime;
t     = 0;

ta = [];           % time
va = [];           % velocities
xa = [];           % position
ka = [];           % kinetic energy
sa = [];           % strain energy

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,1);  % nodal momentum vector
niforce   = zeros(nodeCount,1);  % nodal internal force vector
neforce   = zeros(nodeCount,1);  % nodal external force vector

nvelo     = zeros(nodeCount,1);   % nodal velocity
ndis      = zeros(nodeCount,1);   % nodal displacement
nacce     = zeros(nodeCount,1);   % nodal acceleration

ndis0     = zeros(nodeCount,1);   % old nodal displacement
nvelo0    = zeros(nodeCount,1);   % old nodal velocity
nacce0    = zeros(nodeCount,1);   % old nodal acceleration

activeNodes = unique(elements(unique(pElems),:));
activeDofs  = activeNodes;
K           = zeros(nodeCount,nodeCount); % stiffness matrix
M           = zeros(nodeCount,nodeCount); % mass matrix


%% solution phase
while ( t < time )
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;    
    K(:,:)       = 0; 
    M(:,:)       = 0;
    nacce0(:)    = 0;    
    
    % loop over computational cells or elements
    for e=1:elemCount
        esctr = elements(e,:);
        enode = nodes(esctr);
        mpts  = mpoints{e};        
        % loop over particles
        for p=1:length(mpts)
            pid    = mpts(p);
            vol    = Vp    (pid  );
            m      = Mp    (pid  );
            vel    = vp    (pid,:);
            d      = dis0  (pid,:);
            acc    = acce0 (pid,:);
            xpa    = xp    (pid,:);            
            % shape functions and first derivatives
            N1  = 1 - abs(xpa-enode(1))/deltax;
            N2  = 1 - abs(xpa-enode(2))/deltax;            
            dN1 = -1/deltax;
            dN2 =  1/deltax;            
            % particle mass and momentum to node            
            nmass(esctr(1))      = nmass(esctr(1))     + N1*m;
            nmass(esctr(2))      = nmass(esctr(2))     + N2*m;
            nmomentum(esctr(1))  = nmomentum(esctr(1)) + N1*m*vel;
            nmomentum(esctr(2))  = nmomentum(esctr(2)) + N2*m*vel;  
            nacce0(esctr(1))     = nacce0(esctr(1))    + N1*m*acc;
            nacce0(esctr(2))     = nacce0(esctr(2))    + N2*m*acc;
            % internal force            
            %niforce(esctr(1)) = niforce(esctr(1)) - vol*s(pid)*dN1;
            %niforce(esctr(2)) = niforce(esctr(2)) - vol*s(pid)*dN2;
            % stiffness matrix and mass matrix
            K(esctr(1),esctr(1)) = K(esctr(1),esctr(1)) + vol*dN1*E*dN1;
            K(esctr(1),esctr(2)) = K(esctr(1),esctr(2)) + vol*dN1*E*dN2;
            K(esctr(2),esctr(1)) = K(esctr(2),esctr(1)) + vol*dN2*E*dN1;
            K(esctr(2),esctr(2)) = K(esctr(2),esctr(2)) + vol*dN2*E*dN2;
            
            M(esctr(1),esctr(1)) = M(esctr(1),esctr(1)) + m*N1;
            %M(esctr(1),esctr(2)) = M(esctr(1),esctr(2)) + m*N1*N2;
            %M(esctr(2),esctr(1)) = M(esctr(2),esctr(1)) + m*N2*N1;
            M(esctr(2),esctr(2)) = M(esctr(2),esctr(2)) + m*N2;
        end
    end
    
    nvelo0(activeNodes) = nmomentum(activeNodes)./nmass(activeNodes);
    nacce0(activeNodes) = nacce0(activeNodes)   ./nmass(activeNodes);
    %nvelo0(1)=0;nacce0(1)=0;
    
    dtime2    = dtime*dtime;
    fac       = 1/(beta*dtime2);    
    dtil      = dtime*nvelo0 + 0.5*dtime2*(1-2*beta)*nacce0;    
    A         = fac*M + K;    
    f         = fac*M*dtil; 
    
    % apply Dirichlet BCs
%     bwt       = mean(diag(A));
%     f(1)      = 0;
%     A(1,:)    = 0;
%     A(:,1)    = 0;
%     A(1,1)    = bwt;
    
    ndis     = A\f; 
    
    % update nodal acceleration and velocity    
    nacce = fac*(ndis-dtil);
    nvelo = nvelo0 + (1-gamma)*dtime*nacce0 + gamma*dtime*nacce; 
    ndis(1)=0;
    nvelo(1) = 0; nacce(1)=0;
    
    % update particle velocity and position and stresses
   
    k = 0; u = 0;
    for e=1:elemCount
        esctr = elements(e,:);
        enode = nodes(esctr);
        mpts  = mpoints{e};        
        % loop over particles        
        for p=1:length(mpts)
            pid  = mpts(p);
            
            N1  = 1 - abs(xp(pid)-enode(1))/deltax;
            N2  = 1 - abs(xp(pid)-enode(2))/deltax;            
            dN1 = -1/deltax;
            dN2 =  1/deltax;
            
            d   =  N1*ndis(esctr(1))  + N2*ndis(esctr(2));
            a   =  N1*nacce(esctr(1)) + N2*nacce(esctr(2));
            
            vp(pid)   = vp(pid) +  gamma*dtime*a + (1-gamma)*dtime*acce0(pid);
            xp(pid)   = xp(pid) +  d;
            acce(pid) = a;
            dis(pid)  = d;
            
            Lp          = dN1 * nvelo(esctr(1)) + dN2 * nvelo(esctr(2));
            
            Fp(pid) = (1 + Lp*dtime)*Fp(pid);
            Vp(pid) = Fp(pid)*Vp0(pid);
            dEps    = dtime * Lp;
            s(pid)  = s(pid)   + E * dEps;
            eps(pid)= eps(pid) + dEps;
            
            k    = k + 0.5*vp(pid)^2*Mp(pid);   
            u    = u + 0.5*Vp(pid)*s(pid)*eps(pid);
        end
    end
   
    
    dis0  = dis;
    acce0 = acce;
    
    %nacce0 = nacce;
    
%     for e=1:elemCount
%         esctr = elements(e,:);
%         enode = nodes(esctr,:);
%         mpts  = mpoints{e};        
%         for p=1:length(mpts) % loop over particles
%             pid  = mpts(p);              
%             k    = k + 0.5*vp(pid)^2*Mp(pid);   
%             u    = u + 0.5*Vp(pid)*s(pid)*eps(pid);
%         end
%     end

    % store time,velocty for plotting
    
    cv = 1/sum(Mp)*(dot(Mp,vp));
    ta = [ta;t];    
    va = [va;cv];
    ka = [ka;k];
    sa = [sa;u];
    
    % update the element particle list
    
    pe = floor(xp/deltax)+1;
    
    for e=1:elemCount
        id  = find(pe==e);
        mpoints{e}=id;
    end
    
    activeNodes = unique(elements(unique(pElems),:));
    activeDofs  = activeNodes;

    % advance to the next time step
    
    t = t + dtime;
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

