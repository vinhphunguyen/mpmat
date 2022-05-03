% This file implements the MPM.
%
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% USL with explicit time integration.
%
% Stress wave propagation in a bar. Material models:
% Two scale enhanced damage model of Giang Dinh Nguyen, EFM, 2014.
%
% Use mex functions to speed up the calculations.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% July 2014.

%%

addpath ../fem_util/
addpath ../fem-functions/
addpath ../post-processing/
addpath mex/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties

E   = 20e3;      % Young's modulus
nu  = 0.0;       % Poisson ratio
rho = 2e-8;      % density
ft  = 3.;        % tensile strength
ki  = 1.5e-4;    % damage threshold
alpha=0.95;      % constant in damage evolution law
beta = 4000;     % constant in damage evolution law
h   = 0.5;       % internal length scale
fb  = -0.6*ft;   % distributed force N/mm

stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
De          = elasticityMatrix(E,nu,stressState);

matProp.young = E;
matProp.De    = De;
matProp.ki    = ki;
matProp.alpha = alpha;
matProp.beta  = beta;

tic;

%% Computational grid

ghostCell= 0;
lx       = 100;
ly       = 10;
numx2    = 100;       % number of elements along X direction
numy2    = 1;         % number of elements along Y direction
[mesh]   = buildGrid2D(lx,ly,numx2,numy2, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;
H         = mesh.deltax;

eta=h/H;
%% generate material points

ppc    = 2; % # of particle per cell is ppc x ppc
useGPs = 1; % particles at Gauss points
useGPs = 1; % particles at geometric points regularly divide the elements

monitoredP = 51;

if (useGPs==0)
  dx = mesh.deltax/(ppc);
  dy = mesh.deltay/(ppc);
end

[W,Q]=quadrature(  ppc, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

volume = []; mass   = []; coord  = [];

for e=1:elemCount                 % start of element loop
  sctr = element(e,:);          %  element scatter vector
  pts  = node(sctr,:);
  if (useGPs)
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
  else
    x1 = pts(1,:); % first corner of the cell
    for i=1:ppc
      for j=1:ppc
        x(1) = x1(1) + dx*0.5 + (j-1)*dx;
        x(2) = x1(2) + dy*0.5 + (i-1)*dy;
        volume  = [volume;dx*dy];
        mass    = [mass; dx*dy*rho];
        coord   = [coord;x];
      end
    end
  end
end

pCount = length(volume);

% stored in body1 structure

body1.volume  = volume;
body1.volume0 = volume;
body1.mass    = mass;
body1.coord   = coord;
% if use body1.coord0=coord, then it seems that in MEX function,
% changes made to body1.coord are also present in coord0!!!
coord0  = coord;  % initial particle positions stored to compute u
body1.deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress  = zeros(pCount,3);                % stress
body1.strain  = zeros(pCount,3);                % strain
body1.strain0 = zeros(pCount,3);                % strain0
body1.strainI = zeros(pCount,3);                % strain of localisation band
body1.velo    = zeros(pCount,2);               % velocity
body1.C       = De;
body1.kappa   = zeros(pCount,1);               % max. eqv. strain
body1.damage  = zeros(pCount,1);               % damage variables
body1.status  = zeros(pCount,1);               % localised or not yet
body1.damAtLoc= zeros(pCount,1);               % damage at moment of localisation
%body1.matProp = matProp;

bodies    = cell(1);
bodies{1} = body1;
bodyCount = length(bodies);


%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

bodies = findActiveElemsAndNodes(bodies,mesh);

%% node quantities

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,2);  % nodal momentum vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector
nvelo     = zeros(nodeCount,2);  % nodal velocities
nacce = zeros(nodeCount,2);      % nodal acceleration

%% find boundary conditions

eps=1e-12;
left   = find(abs(node(:,1))   <eps);
right  = find(abs(node(:,1)-lx)<eps);

ff = ly/numy2;
neforce(right([1 end]),1) = ff*fb/2.;
neforce(right([2:end-1]),1) = ff*fb;

fixedBoth = left;


xxx = coord(1,2);
row = find(abs(coord(:,2)-xxx)<eps);

% indexes of particles whose quantities are monitored for post-processing
% this is dependent on the discretization!!!

mid = find(abs(coord(:,1)-50.2113)<1e-4);
rig = find(abs(coord(:,1)-99.7887)<1e-4);

lef = 1;
mid = mid(1);
rig = rig(1);

%% plot mesh, particles

figure
set(gca,'FontSize',14)
hold on
plot_mesh(node,element,'Q4','k-',1.6);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','cy-',1.6);
xp1 = bodies{1}.coord;
plot(xp1(:,1),xp1(:,2),'k.','markersize',11);
plot(xp1(row,1),xp1(row,2),'r.','markersize',15);
plot(node(left,1),node(left,2),'cy.','markersize',15);
plot(node(right,1),node(right,2),'cy.','markersize',15);
%plot(node(bodies{1}.nodes,1),node(bodies{1}.nodes,2),'b.','markersize',15);
%axis off

ta = [0];           % time
sa = [0 0 0];           % recorded stress
va = [0 0 0];

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

cfl   = 0.05;
c     = sqrt(E/rho);
dtime = (mesh.deltax/c)*cfl;
time  = 2.*lx/c;
ftime = time/16;


t     = 0;


nsteps = floor(time/dtime);

str     = cell(nsteps,1);
eqv     = cell(nsteps,1);
istep = 1;

nx       = 1; ny = 0; % hard coded
n        = zeros(3,2);
n(1,1)   = nx;
n(2,2)   = ny;
n(3,1)   = ny;
n(3,2)   = nx;

while ( t < time )
  disp(['time step ',num2str(t)])
  
  [nmass,nmomentum,niforce] = ParticlesToNodes(bodies,mesh);
  
  % update nodal momenta
  
  if ( t > ftime)
    neforce(right([1 end]),1)   = 0.;
    neforce(right([2:end-1]),1) = 0.;
  end
  
  nforce = neforce + niforce;
  
  activeNodes=bodies{1}.nodes;
  nmomentum(activeNodes,:) = nmomentum(activeNodes,:) + nforce(activeNodes,:)*dtime;
  
  nvelo(activeNodes,1) = nmomentum(activeNodes,1)./nmass(activeNodes);
  nvelo(activeNodes,2) = nmomentum(activeNodes,2)./nmass(activeNodes);
  
  nacce(activeNodes,1) = nforce(activeNodes,1)./nmass(activeNodes);
  nacce(activeNodes,2) = nforce(activeNodes,2)./nmass(activeNodes);
  
  % boundary conditions
  nvelo(fixedBoth,:) = 0; % fixed boundary conditions
  nacce(fixedBoth,:)    = 0; % fixed boundary conditions
  
  %   for ib=1:bodyCount
  %     body      = bodies{ib};
  %     elems     = body.elements;
  %     mpoints   = body.mpoints;
  %     for ie=1:length(elems)         % loop over computational cells or elements
  %       e     = elems(ie);
  %       esctr = element(e,:);    % element connectivity
  %       mpts  = mpoints{e};        % particles inside element e
  %       for p=1:length(mpts)       % loop over particles
  %         pid  = mpts(p);
  %         xp   = body.coord(pid,:);
  %         xp0  = body.coord(pid,:);
  %         vp   = body.velo(pid,:);
  %         Lp   = zeros(2,2);
  %         for i=1:length(esctr)
  %           id  = esctr(i);
  %           x   = xp0 - node(id,:);
  %           [N,dNdx]=getMPM2D(x,mesh.deltax,mesh.deltay);
  %           if (nmass(id)~=0)
  %             xp  = xp + dtime * N*nvelo(id,:);
  %             vp  = vp + dtime * N*nacce(id,:);
  %             Lp  = Lp + nvelo(id,:)'*dNdx;% particle gradient velocity
  %           end
  %         end
  %
  %         bodies{ib}.coord(pid,:)= xp;
  %         bodies{ib}.velo(pid,:) = vp;
  %
  %         F       = ([1 0;0 1] + Lp*dtime)*reshape(bodies{ib}.deform(pid,:),2,2);
  %         bodies{ib}.deform(pid,:) = reshape(F,1,4);
  %         bodies{ib}.volume(pid  ) = det(F)*bodies{ib}.volume0(pid);
  %         dEps    = dtime * 0.5 * (Lp+Lp');
  %         bodies{ib}.strain(pid,:)  = bodies{ib}.strain(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
  %
  %         dsigma  = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
  %         bodies{ib}.stress(pid,:)  = bodies{ib}.stress(pid,:) + dsigma';
  %       end
  %     end
  %   end
  
  bodies{1}.strain0 = bodies{1}.strain;
  
  UpdateParticles(bodies,mesh,nvelo,nacce,dtime); % MEX function
  
  
  for ib=1:bodyCount
    body = bodies{ib};
    for p=1:length(bodies{ib}.mass) % loop over computational cells or elements
      kappa0  = body.kappa(p,:);    % kappa at time t 
      status  = body.status(p);     % status at time t 
      strain  = body.strain(p,:)';  % total strain at time t+dt
      strain0 = body.strain0(p,:)'; % total strain at time t
      if (status==0)                % outside of localisation band
        eqvS = sqrt((1./E)*strain'*De*strain);
        f   = eqvS - kappa0;
        
        kappa = kappa0;
        
        if ( f >= 0.)
          kappa = eqvS;
        end
        damage = 0.;
        if ( kappa >= ki )
          damage = 1. - (ki/kappa)*(1.-alpha+alpha*exp(-beta*(kappa-ki)));
        end
        
        stress = (1.-damage)*De*strain;
        
        % check localisation occurence
        if (damage>0)
          bodies{ib}.status(p)   = 1;
          bodies{ib}.damAtLoc(p) = damage;
          bodies{ib}.strainI(p,:)= strain;
        end
      else
        % compute ao
        damL = bodies{ib}.damAtLoc(p);
        ao   = (1-damL)*De;
        % compute ai
        straini = bodies{ib}.strainI(p,:)';
        eqvS    = sqrt((1./E)*straini'*De*straini);
        f       = eqvS - kappa0;
        
        kappa = kappa0;
        
        if ( f >= 0.)
          kappa = eqvS;
        end
        
        damage = 0.;
        if ( kappa >= ki )
          damage = 1. - (ki/kappa)*(1.-alpha+alpha*exp(-beta*(kappa-ki)));
        end
  
        a        = ki / (kappa * kappa);
        b        = a * kappa;
        c        = alpha * exp(-beta * (kappa-ki) );
        dddk     = a * (1.0 - alpha+ c) + beta * b * c;
        
        deqvdeps = (1/E/eqvS)*De*straini;
        ai       = (1-damage)*De - dddk*De*straini*deqvdeps';
        % acoustic tensors
        Ao   = n'*ao*n;
        Ai   = n'*ai*n;
        C        = (eta/h)*Ao + (1-eta)/h*Ai;
        dstrain  = (strain - strain0);
        djumpu   = inv(C)*n'*ao*dstrain;
        depsi    = (1/h)*n*djumpu;
        % update strain at the localisation band
        bodies{ib}.strainI(p,:) = bodies{ib}.strainI(p,:) + depsi';
        % compute stress increment at outside (also the total)
        depso    = dstrain - eta*depsi;
        dsigma   = (1/(1-eta))*ao*depso;
        stress   = bodies{ib}.stress(p,:)' + dsigma;  
      end
      
      bodies{ib}.stress(p,:) = stress';
      bodies{ib}.kappa(p,:)  = kappa;
    end
  end
  
  
  bodies = findActiveElemsAndNodes(bodies,mesh);
  
  str{istep} = bodies{1}.stress;
  eqv{istep} = bodies{1}.kappa;
  
  % advance to the next time step
  t = t + dtime;
  istep = istep + 1;
  
  s  = bodies{1}.stress(monitoredP,1);
  u  = bodies{1}.coord(pCount,1)-coord0(pCount,1);
  ta = [ta;t];
  sa = [sa;bodies{1}.stress(lef,1) bodies{1}.stress(mid,1) bodies{1}.stress(rig,1)];
  va = [va;bodies{1}.velo(lef,1) bodies{1}.velo(mid,1) bodies{1}.velo(rig,1)];
end

%%

sig =str{900};
sig900=sig(row,1);
sig =str{500};
sig500=sig(row,1);
sig =str{200};
sig200=sig(row,1);

row=[];

for i=1:elemCount/5
  e   = 1 + 4*(i-1);
  row = [row; e; e+2];
end


kap =eqv{1900};
kap2000=kap(row,1);
kap =eqv{3000};
kap3000=kap(row,1);
kap =eqv{4000};
kap3100=kap(row,1);
xx  = coord0(row,1);
%
figure
set(gca,'FontSize',14)
hold on
plot(xx,kap2000,'red-','LineWidth',1.1);
plot(xx,kap3000,'blue*','LineWidth',1.1);
plot(xx,kap3100,'black--','LineWidth',1.1);
xlabel('distance')
ylabel('stress');
legend('200','500','900');
%%

%%
figure
set(gca,'FontSize',14)
hold on
plot(ta*1e6,sa(:,1),'r-','LineWidth',1.1);
plot(ta*1e6,sa(:,2),'b-','LineWidth',1.1);
plot(ta*1e6,sa(:,3),'black-','LineWidth',1.1);
xlabel('time [microsecond]')
ylabel('stress');
legend('x=0', 'x=0.5L', 'x=L');
set(gca,'XTick',[0:40:400])
xlim([0, 400])

%%

figure
set(gca,'FontSize',14)
hold on
plot(ta*1e6,va(:,1),'r-','LineWidth',1.1);
plot(ta*1e6,va(:,2),'b-','LineWidth',1.1);
plot(ta*1e6,va(:,3),'black-','LineWidth',1.1);
xlabel('time [microsecond]')
ylabel('stress');
legend('x=0', 'x=0.5L', 'x=L');
set(gca,'XTick',[0:40:400])
xlim([0, 400])

disp([num2str(toc),'   DONE '])