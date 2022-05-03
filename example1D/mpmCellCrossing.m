% This file implements the Material Point Method of Sulsky 1994.
% One dimensional problem with material points.
% The grid is one two-noded linear element.
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% August 2015.

%%

addpath ../grid/
addpath ../basis/
addpath ../particleGen/
addpath ../constitutiveModels/
addpath ../util/
addpath ../geoMesh/
addpath ../postProcessing/


%%
clc
clear all

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%%
%
E      = 1000;               % Young modulus
rho    = 1000;               % density
c      = sqrt(E/rho);
v      = 0.1;

%%
%  Computational grid: two-noded elements
L     = 1;
ne    = 5;
[mesh]=buildGrid1D(L,ne,0);

elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
elements  = mesh.element;
nodes     = mesh.node;
deltax    = mesh.deltax;
omegaC    = deltax;

%%
% Material points:
ppc       = 1; % # particles per cell
[mesh1]   = buildGrid1D(deltax,1,0);
particles = buildParticlesGeom(mesh1,ppc,rho);

pCount = 2*particles.pCount;
xp   = [particles.xp;L-particles.xp];
xp0  = xp;
vp   = [particles.vp;particles.vp];
Vp   = [particles.Vp;particles.Vp];
Vp0  = Vp;
Fp   = ones(pCount,1);
s    = [particles.s;particles.s];
eps  = [particles.eps;particles.eps];
Mp   = [particles.Mp;particles.Mp];

% initial velocities, initial stress=0
for p=1:pCount
  if xp(p) < L/2
    vp(p) = v;
  else
    vp(p) = -v;
  end
end
vp0=vp;

%%
figure
hold on
plot(nodes,zeros(nodeCount,1)+1/2,'r-s');
plot(xp,zeros(pCount,1)+1/2,'b*');
%axis([0 L+0.5 0 1.])

%%
% data structure to store the material points for each element
% this data structure is updated for every time step

pElems  = ones(pCount ,1);
mpoints = cell (elemCount ,1);

for p=1:pCount
    x = xp(p);
    e = floor(x/deltax) + 1; % MPM
    
    pElems(p) = e; % particle "p" stays in element "e"
    for e=1:elemCount
        id = find(pElems==e);
        mpoints{e}=id ; % mpoints{e}?> indices of particles in "e"
    end
end

%% Time loop

tol = 1e-8;

dtime = 0.001;
time  = 4;
t     = 0;

ta = [];           % time
va = [];           % velocities
xa = [];           % position
ka = [];           % kinetic energy
sa = [];           % strain energy
sta= [];           % stress 

force=[];
strain=[];

pos={};
ix     = 1;
istep = 1;

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,1);  % nodal momentum vector
niforce   = zeros(nodeCount,1);  % nodal internal force vector
neforce   = zeros(nodeCount,1);  % nodal external force vector

while ( t < time )
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;    
    % loop over computational cells or elements
    for e=1:elemCount
        esctr = elements(e,:);
        enode = nodes(esctr);
        mpts  = mpoints{e};        
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            % shape functions and first derivatives
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
    nmomentum = nmomentum + niforce*dtime;
    
    % update particle velocity and position and stresses
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
            
            if nmass(esctr(1)) > tol
                vp(pid)  = vp(pid) + dtime * N1*niforce(esctr(1))/nmass(esctr(1));
            end
            
            if nmass(esctr(2)) > tol
                vp(pid)  = vp(pid) + dtime * N2*niforce(esctr(2))/nmass(esctr(2));
            end
            
            xp(pid)  = xp(pid) + dtime * (N1*nmomentum(esctr(1))/nmass(esctr(1)) + ...
                N2*nmomentum(esctr(2))/nmass(esctr(2)));
            
            v1 = nmomentum(esctr(1))/nmass(esctr(1));
            v2 = nmomentum(esctr(2))/nmass(esctr(2));
            
            %if ( esctr(1) == 1 ) v1 = 0; end
            % gradient velocity
            
            Lp = dN1 * v1 + dN2 * v2;
            
            Fp(pid) = (1 + Lp*dtime)*Fp(pid);
            Vp(pid) = Fp(pid)*Vp0(pid);
            dEps    = dtime * Lp;
            s(pid)  = s(pid)   + E * dEps;
            eps(pid)= eps(pid) + dEps;
            
            k = k + 0.5*vp(pid)^2*Mp(pid);
            u = u + 0.5*s(pid)*eps(pid)*Vp(pid);
        end
    end
    
        if ( t > 2.55 ) 
      force=[force niforce]; 
      strain=[strain eps']; 
    end
    
%     if (xp(1)>2*deltax)
%       disp('passing xxxx')   
%       nmomentum.\nmass
%       nmass
%       niforce
%       %pause
%     end

    
    % store time,velocty for plotting
        
    ta = [ta;t];        
    ka = [ka;k];
    sa = [sa;u];
    sta = [sta;s(1)];
    
    if  ( mod(istep,100) == 0 )
    pos{ix} = xp; ix = ix+1;
    end
    
    % update the element particle list
    
    pe = floor(xp/deltax)+1;
    
    
    for e=1:elemCount
        id  = find(pe==e);
        mpoints{e}=id;
    end

    % advance to the next time step
    
    t = t + dtime;
    istep = istep + 1;
end

%%

% exact solution

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
%ka = ka*1000;
%sa = sa*1000;
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
%axis([2.53 2.57 0 3])

%savefile = 'mpmCollision1D.mat';
%save(savefile,'ta','ka','sa');
%%
for i=1:length(pos)
  figure
  hold on
  plot(nodes,zeros(nodeCount,1)+1/2,'r-s');
  plot(pos{i},zeros(pCount,1)+1/2,'b*');
  title( sprintf('%i',num2str(90*dtime*i) ));
  frame(i) = getframe;
end
close all;
movie(frame);

