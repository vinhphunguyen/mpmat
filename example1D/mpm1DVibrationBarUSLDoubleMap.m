% This file implements the Material Point Method of Sulsky 1994.
% One dimensional problem with material points.
% The grid is one two-noded linear element.
% This example is taken from "Caveats on the Implementation of the Generalized
% Material Point Method", Buzzi et al, CMES 2008.
%
% Update stress last and double mapping technique.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% June 2013.

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
n      = 10;            % mode number
c      = sqrt(E/rho);
beta1  = (2*n-1)*0.5*(pi/L);
omega1 = beta1*c;

%%
%  Computational grid: two-noded elements

[mesh]=buildGrid1D(L,100,0);

elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
elements  = mesh.element;
nodes     = mesh.node;
deltax    = mesh.deltax;

%%
% Material points:
ppc = 2; % # particles per cell
particles=buildParticles(mesh,ppc,rho);

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


dtime = 0.01*deltax/c;
time  = (2*pi/omega1)*10;
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
nvelo     = zeros(nodeCount,1);  % nodal external force vector

while ( t < time )
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    nvelo(:)     = 0;
    % loop over computational cells or elements
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
            
            % particle mass and momentum to node
            
            nmass(esctr(1))      = nmass(esctr(1))     + N1*Mp(pid);
            nmass(esctr(2))      = nmass(esctr(2))     + N2*Mp(pid);
            nmomentum(esctr(1))  = nmomentum(esctr(1)) + N1*Mp(pid)*vp(pid);
            nmomentum(esctr(2))  = nmomentum(esctr(2)) + N2*Mp(pid)*vp(pid);
            
            % internal force
            
            niforce(esctr(1)) = niforce(esctr(1)) - Vp(pid)*s(pid)*dN1;
            niforce(esctr(2)) = niforce(esctr(2)) - Vp(pid)*s(pid)*dN2;
        end
    end
    
    % update nodal momenta
    
    nmomentum(1)  = 0; % Boundary conditions f1 = m1*a1, a1=0
    niforce(1)    = 0;
        
    nmomentum = nmomentum + niforce*dtime;
        
    % update particle velocity, double mapping to get nodal velocities
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
                        
            vp(pid)  = vp(pid) + dtime * N1*niforce(esctr(1))/nmass(esctr(1));                
            vp(pid)  = vp(pid) + dtime * N2*niforce(esctr(2))/nmass(esctr(2));                

            nvelo(esctr(1))  = nvelo(esctr(1))  + N1*Mp(pid)*vp(pid);
            nvelo(esctr(2))  = nvelo(esctr(2))  + N2*Mp(pid)*vp(pid);
        end
    end
    
    nvelo = nvelo./nmass;
    
    k = 0;
    u = 0;
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
            % updated positions                                    
            xp(pid)  = xp(pid) + dtime * (N1*nmomentum(esctr(1))/nmass(esctr(1)) + ...
                 N2*nmomentum(esctr(2))/nmass(esctr(2)));

            v1 = nvelo(esctr(1));
            v2 = nvelo(esctr(2));
    
            if ( esctr(1) == 1 ) v1 = 0; end
            % gradient velocity
            
            Lp = dN1 * v1 + dN2 * v2;
            
            Fp(pid) = (1 + Lp*dtime)*Fp(pid);
            Vp(pid) = Fp(pid)*Vp0(pid);                                    
            dEps    = dtime * Lp;                                    
            s(pid)  = s(pid)   + E * dEps;    
            eps(pid)= eps(pid) + dEps;  
            
            k = k + 0.5*vp(pid)^2*Mp(pid);
            u = u + 0.5*s(pid)*eps(pid)*Mp(pid);
        end
    end
           
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
    
    % advance to the next time step
    
    t = t + dtime;
end

%%

% exact solution

%vExact = 0.1*cos(omega1.*ta)*sin(beta1*0.5*L);
vExact = 0.1/(beta1*L)*cos(omega1.*ta);
uExact = (0.1/omega1)*sin(omega1.*ta)*sin(beta1*0.5*L);

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
% %axis([0 100 -0.15 0.2])

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

