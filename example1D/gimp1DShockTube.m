% This file implements the Generalized Material Point Method.
% One dimensional problem with material points.
% The grid is one two-noded linear element.
%
% Sod's shock tube problem in gas dynamics.
%
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia.
% June 2015.

%%

addpath ../fem_util/


%%
clc
clear all

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

% Material properties

gamma = 1.4;   % perfect gas
c0    = 2.;
c1    = 1.;
viscosity = 1; % switch on/off artifical viscosity

%%
%  Computational grid: two-noded elements
elemCount = 200;
L = 1;
[mesh]=buildGrid1D(L,elemCount,1);

elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
elements  = mesh.element;
nodes     = mesh.node;
deltax    = mesh.deltax;

%%
% Material points:
ppc = 3; % # particles per cell

Mp   = [];
Vp   = [];
xp   = [];
sp   = [];
rp   = [];   % density
ep   = [];   % internal energy

dx   = mesh.deltax/(ppc+1);
vol  = mesh.deltax/(ppc);

rho  = 1;
p    = 1;

for e=1:mesh.elemCount                 % start of element loop
    sctr = mesh.element(e,:);          %  element scatter vector
    pts  = nodes(sctr,:);
    % skip first and last element (ghost cells)
    if (e==1) || (e==mesh.elemCount)
       continue;
    end
    for q=1:ppc                           
        x   = pts(1) + vol*0.5 + (q-1)*vol; 
        if ( x > 0.5 )
          rho = 0.125; 
          p   = 0.1;
        end
        Vp  = [Vp;vol];
        Mp  = [Mp;vol*rho];
        xp  = [xp;x];
        sp  = [sp;-p];       % stress = -pressure;
        rp  = [rp;rho];
        ep  = [ep;p/(gamma-1)/rho];
    end
end

pCount = length(xp);
Vp0 = Vp;
Fp  = ones(pCount,1);
eps = zeros(pCount,1);
vp  = zeros(pCount,1);               % velocity

%%
hold on
plot(nodes,zeros(nodeCount,1)+1/4,'r-s','MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',9,'LineWidth',1.1);
plot(xp,zeros(pCount,1)+1/4,'bo','MarkerEdgeColor','k',...
                    'MarkerFaceColor','r',...
                    'MarkerSize',9,'LineWidth',1.1);
% axis([-10 35 0 1/2])


%%
% data structure to store the material points for each element
% this data structure is updated for every time step

pElems  = ones(pCount ,1);
mpoints = cell (elemCount ,1);

for p=1:pCount
    x = xp(p);
    e = floor(x/deltax) + 2; % + 2 because x starts -1
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

dtime = 0.0001;
time  = 0.143;
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
    disp(['time step ',num2str(t)])
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    nvelo(:)     = 0;
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
                [phi,dphi]=getGIMP(x,deltax,lp);
                % particle mass and momentum to node
                nmass(id)     = nmass(id)     + phi*Mp(pid);
                nmomentum(id) = nmomentum(id) + phi*Mp(pid)*vp(pid);
                niforce(id)   = niforce(id)   - Vp(pid)*sp(pid)*dphi;
            end
        end
    end
    
    % update nodal momenta
    
    % be careful with left/right ghost cells
    
    niforce(1)    = 0; niforce(2)    = 0;    
    nmomentum(1)  = 0; nmomentum(2)  = 0; 
    niforce(nodeCount)    = 0; niforce(nodeCount-1)    = 0;
    nmomentum(nodeCount)  = 0; nmomentum(nodeCount-1)  = 0; 
    
    nmomentum = nmomentum + niforce*dtime;
 
    % update particle velocity and position and stresses
   
    for e=1:elemCount
        esctr = gimpElement{e};
        mpts  = mpoints{e};
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
                   
            for i=1:length(esctr)
                id    = esctr(i);
                x     = xp(pid,:) - nodes(id,:);
                % shape functions and first derivatives
                [phi,dphi]=getGIMP(x,deltax,lp);
                
                if nmass(id) > 0
                    vp(pid) = vp(pid) + dtime * phi*niforce(id)/nmass(id);                                                                                
                end                                                                                                
            end     
            
            for i=1:length(esctr)
              id    = esctr(i);
              x     = xp(pid,:) - nodes(id,:);
              % shape functions and first derivatives
              [phi,dphi]=getGIMP(x,deltax,lp);                            
              nvelo(id)  = nvelo(id)  + phi*Mp(pid)*vp(pid);            
            end
        end
    end
    
    for i=1:length(nvelo)
      if nmass(i) > 0
        nvelo(i) = nvelo(i)/nmass(i);
      end
    end
    
    nvelo(1)          = 0; nvelo(2)          = 0;
    nvelo(nodeCount)  = 0; nvelo(nodeCount-1)= 0;
    
    for e=1:elemCount
        esctr = gimpElement{e};
        mpts  = mpoints{e};
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            Lp = 0;            
            for i=1:length(esctr)
                id    = esctr(i);
                x     = xp(pid,:) - nodes(id,:);
                % shape functions and first derivatives
                [phi,dphi]=getGIMP(x,deltax,lp);   
                
                if nmass(id) > 0                                                           
                    xp(pid) = xp(pid) + dtime * phi*nmomentum(id)/nmass(id);                                        
                end                                                                
                % gradient velocity                
                Lp = Lp + dphi * nvelo(id);
            end
            Fp(pid) = (1 + Lp*dtime)*Fp(pid);
            Vp(pid) = Fp(pid)*Vp0(pid);
            
            % update particle density using continuity equation
            rp(pid) = rp(pid)/(1+dtime*Lp);
            % update internal energy
            ee0     = ep(pid) + dtime*sp(pid)*Lp/rp(pid);
            % compute pressure using EOS
            pres = ( (gamma-1)*rp(pid)*ep(pid) ) / ( 1 - (gamma-1)*dtime*Lp);
            q    = 0;
            if ( Lp < 0 )
              a = sqrt(gamma*pres/rp(pid));
              q = rp(pid)*deltax*(c0*deltax*Lp*Lp-c1*a*Lp);
            end
            pres = pres + viscosity*q;
            sig  = -pres;
            ee   = ep(pid) + dtime*sig*Lp/rp(pid);
            ep(pid) = ee;
            sp(pid) = sig;                        
        end
    end
    
    % store time,velocity for plotting
    
  
    
    % update the element particle list
    
    pe = floor(xp/deltax)+2;
    
    for e=1:elemCount
        id  = find(pe==e);
        mpoints{e}=id;
    end
    
    % advance to the next time step
    
    t = t + dtime;
end

%%

% exact solution

exact = load('sod-exact.mat');

figure
subplot(2,2,1)
set(gca,'FontSize',14)
hold on
plot(xp,rp,'b-','LineWidth',1.6);
plot(exact.data.x,exact.data.rho,'r-','LineWidth',1.6);
xlabel('Position')
ylabel('Density')
axis([0 1 0 1.2])
ax = gca;
ax.XTick = [0.0 0.2 0.4 0.6 0.8 1.0];
legend('MPM','Exact')

subplot(2,2,2)
set(gca,'FontSize',14)
hold on
plot(xp,-sp,'b-','LineWidth',1.6);
plot(exact.data.x,exact.data.P,'r-','LineWidth',1.6);
xlabel('Position')
ylabel('pressure')
legend('MPM','Exact')
ax = gca;
ax.XTick = [0.0 0.2 0.4 0.6 0.8 1.0];
axis([0 1 0 1.2])

subplot(2,2,3)
set(gca,'FontSize',14)
hold on
plot(xp,vp,'b-','LineWidth',1.6);
plot(exact.data.x,exact.data.u,'r-','LineWidth',1.6);
xlabel('Position')
ylabel('Velocity')
ax = gca;
ax.XTick = [0.0 0.2 0.4 0.6 0.8 1.0];
axis([0 1 0 1.6])


subplot(2,2,4)
set(gca,'FontSize',14)
hold on
plot(xp,ep,'b-','LineWidth',1.6);
plot(exact.data.x,exact.data.e,'r-','LineWidth',1.6);
xlabel('Position')
ylabel('internal energy')
ax = gca;
ax.XTick = [0.0 0.2 0.4 0.6 0.8 1.0];
axis([0 1 1 3.5])


