% This file implements the Material Point Method of Sulsky 1994.
% One dimensional problem with material points.
% The grid is one two-noded linear element.
%
% Sod's shock tube problem in gas dynamics.
%
% Update stress last and double mapping technique.
%
% Vinh Phu Nguyen
% Monash University, Australia.
% September 2019.

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

% Material properties

gamma = 1.4;   % perfect gas
c0    = 2.;
c1    = 1.;
viscosity = 1; % switch on/off artifical viscosity

%%
%  Computational grid: two-noded elements
L = 1;
[mesh]=buildGrid1D(L,200,0);

elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
elements  = mesh.element;
nodes     = mesh.node; % [0,1] is the domain
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
    e = floor(x/deltax) + 1;
    pElems(p) = e; % particle "p" stays in element "e"
    for e=1:elemCount
        id = find(pElems==e);
        mpoints{e}=id ; % mpoints{e}?> indices of particles in "e"
    end
end


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
nstress   = zeros(nodeCount,1);  % nodal stress


while ( t < time )
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    nvelo(:)     = 0;
    nstress(:)   = 0;
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
            
            niforce(esctr(1)) = niforce(esctr(1)) - Vp(pid)*sp(pid)*dN1;
            niforce(esctr(2)) = niforce(esctr(2)) - Vp(pid)*sp(pid)*dN2;
        end
    end
    
    % update nodal momenta
    
    nmomentum(1)  = 0; nmomentum(end)  = 0;
    niforce(1)    = 0; niforce(end)    = 0;
        
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
    
    nvelo(1) = 0; nvelo(end)=0;
    
    for e=1:elemCount
        esctr = elements(e,:);
        enode = nodes(esctr);        
        mpts  = mpoints{e};        
        % loop over particles  
        doo = 1;
        for p=1:length(mpts)
            pid  = mpts(p);
            
            N1  = 1 - abs(xp(pid)-enode(1))/deltax;
            N2  = 1 - abs(xp(pid)-enode(2))/deltax;
            
            dN1 = -1/deltax;
            dN2 =  1/deltax;
                                                
            xp(pid)  = xp(pid) + dtime * (N1*nmomentum(esctr(1))/nmass(esctr(1)) + ...
                 N2*nmomentum(esctr(2))/nmass(esctr(2)));

            v1 = nvelo(esctr(1));
            v2 = nvelo(esctr(2));
    
            %if ( esctr(1) == 1 ) v1 = 0; end
            % gradient velocity
            
            Lp = dN1 * v1 + dN2 * v2;
            
            Fp(pid) = (1 + Lp*dtime)*Fp(pid);
            Vp(pid) = Fp(pid)*Vp0(pid);                                    
            dEps    = dtime * Lp;    
            
            % update particle density using continuity equation
            rp(pid) = rp(pid)/(1+dtime*Lp);
            % update internal energy
            %ee0     = ep(pid) + dtime*sp(pid)*Lp/rp(pid);
            % compute pressure using EOS
            pres = ( (gamma-1)*rp(pid)*ep(pid) ) / ( 1 - (gamma-1)*dtime*Lp);
            % compute artifical viscosity q
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
               
    for e=1:elemCount
        esctr = elements(e,:);
        enode = nodes(esctr);        
        mpts  = mpoints{e};        
        % loop over particles          
        for p=1:length(mpts)
            pid  = mpts(p);            
            N1  = 1 - abs(xp(pid)-enode(1))/deltax;
            N2  = 1 - abs(xp(pid)-enode(2))/deltax;
            sig = sp(pid);
            nstress(esctr(1)) = nstress(esctr(1)) + N1*Mp(pid)*sig;
            nstress(esctr(2)) = nstress(esctr(2)) + N2*Mp(pid)*sig;
        end
    end
    
    nstress = nstress./nmass;
    
     % store time,velocty for plotting
    
    cv = 1/sum(Mp)*(dot(Mp,vp));
    ta = [ta;t];    
    va = [va;cv];
    
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

%exact = load('sod-exact.mat');

figure
subplot(2,2,1)
set(gca,'FontSize',14)
hold on
plot(xp,rp,'b-','LineWidth',1.6);
%plot(exact.data.x,exact.data.rho,'r-','LineWidth',1.6);
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

%%
figure(2)
subplot(2,2,1)
set(gca,'FontSize',14)
hold on
plot(nodes(2:end-1),nvelo(2:end-1),'b-','LineWidth',1.6);
plot(exact.data.x,exact.data.u,'r-','LineWidth',1.6);
xlabel('Position')
ylabel('Velocity')
ax = gca;
ax.XTick = [0.0 0.2 0.4 0.6 0.8 1.0];
axis([0 1 0 1.6])
title('At grid nodes')
legend('MPM','Exact')


subplot(2,2,2)
set(gca,'FontSize',14)
hold on
plot(xp,vp,'b-','LineWidth',1.6);
plot(exact.data.x,exact.data.u,'r-','LineWidth',1.6);
xlabel('Position')
ylabel('Velocity')
ax = gca;
ax.XTick = [0.0 0.2 0.4 0.6 0.8 1.0];
axis([0 1 0 1.6])
title('At material points')

figure(2)
subplot(2,2,4)
set(gca,'FontSize',14)
hold on
plot(xp,-sp,'b-','LineWidth',1.6);
plot(exact.data.x,exact.data.P,'r-','LineWidth',1.6);
ylabel('Pressure')
xlabel('Position')
ax = gca;
ax.XTick = [0.0 0.2 0.4 0.6 0.8 1.0];
axis([0 1 0 1.6])
title('At material points')

subplot(2,2,3)
set(gca,'FontSize',14)
hold on
plot(nodes(2:end-1),-nstress(2:end-1),'b-','LineWidth',1.6);
plot(exact.data.x,exact.data.P,'r-','LineWidth',1.6);
ylabel('Pressure')
xlabel('Position')
ax = gca;
ax.XTick = [0.0 0.2 0.4 0.6 0.8 1.0];
axis([0 1 0 1.6])
title('At  grid nodes')



