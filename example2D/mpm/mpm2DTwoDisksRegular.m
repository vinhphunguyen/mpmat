% This file implements the Material Point Method of Sulsky 1994.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% Regular initial particle distribution. 
%
% Shape functions: using standard FE shape function routine by converting
% particle position to natural coordinates. 
%
% Two elastic disks come into contact.
% Example taken from Sulsky et al. 1994 paper.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% June 2013.
% August 2015: map particles' stress to grid 

%%

addpath ../../grid/
addpath ../../basis/
addpath ../../particleGen/
addpath ../../constitutiveModels/
addpath ../../util/
addpath ../../geoMesh/
addpath ../../externals/PolyMesher/
addpath ../../externals/
addpath ../../postProcessing/

%%
clc
clear 
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties 
% 

E   = 1000;        % Young's modulus
nu  = 0.3;         % Poisson ratio
rho = 1000;        % density
kappa = 3-4*nu;    % Kolosov constant
mu    = E/2/(1+nu);% shear modulus

fac   = 1; % 1: explicit time 
             % 100: implicit is better 
v     = 0.1/fac;   % initial particle velocity 

interval     = 100;% time interval for saving vtp files.
vtkFileName  = 'mpm2DTwoDisks';

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])



%% Material properties
E      = 1000;              % Young modulus
nu     = 0.3;               % Poisson ratio
rho    = 1000;              % density
c      = sqrt(E/rho);
v0     = 0.1;
stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);

interval     = 100;% time interval for saving vtp files.
vtkFileName  = 'mpm2DTwoDisks';

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell=0;
lx      = 1;
ly      = 1;
noX0    = 20;      % number of elements along X direction
noY0    = noX0;        % number of elements along Y direction
[mesh]  = buildGrid2D(lx,ly,noX0,noY0, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;
omegaC    = mesh.deltax*mesh.deltay;
deltax    = mesh.deltax;
deltay    = mesh.deltay;

%% generate material points
ppc           = [2 2];
cir1.center   = [0.2 0.2];
cir1.radius   = 0.2;
cir2.center   = [0.8 0.8];
cir2.radius   = 0.2;
[res1]        = generateMPForCircle(cir1,ppc,mesh);
[res2]        = generateMPForCircle(cir2,ppc,mesh);

pCount  = 2*size(res1.position,1);
Vp      = [res1.volume;res2.volume];
Vp0     = Vp;
Mp      = Vp*rho;
xp      = [res1.position;res2.position];
coords0 = xp;
Fp      = repmat([1 0 0 1],pCount,1);     % gradient deformation
s       = zeros(pCount,3);                % stress
eps     = zeros(pCount,3);                % strain
vp      = zeros(pCount,2);                % velocity


% initial velocities, initial stress=0
for p=1:pCount
  if xp(p,1) < 0.5
    vp(p,:) = [v0 v0];
  else
    vp(p,:) = [-v0 -v0];
  end
end

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

pElems  = ones(pCount,1);
mpoints = cell(elemCount,1);

for p=1:pCount
    x = xp(p,1);
    y = xp(p,2);
    e = floor(x/deltax) + 1 + mesh.numx*floor(y/deltay);
    pElems(p) = e;
end

for e=1:elemCount
    id  = find(pElems==e);
    mpoints{e}=id;
end


%% node quantities

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,2);  % nodal momentum vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector (no need for this exam)

%% plot mesh, particles

figure(1)
hold on
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(xp(:,1),xp(:,2),'k.','markersize',10);
axis off

ta = [];           % time
ka = [];           % kinetic energy 
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.001;
%time  = 30*dtime; % low velocity, used to test implicit
time  = 4% explicit, high velocity 
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
    disp(['time step ',num2str(t)])
    % reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;    
    % loop over computational cells or elements
    for e=1:elemCount
        esctr = element(e,:);      % element connectivity
        enode = node(esctr,:);     % element node coords        
        mpts  = mpoints{e};        % particles inside element e         
        % loop over particles        
        for p=1:length(mpts)
            pid  = mpts(p);
            pt(1)= (2*xp(pid,1)-(enode(1,1)+enode(2,1)))/deltax;
            pt(2)= (2*xp(pid,2)-(enode(2,2)+enode(3,2)))/deltay;
            
            [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
            J0       = enode'*dNdxi;             % element Jacobian matrix
            invJ0    = inv(J0);
            dNdx     = dNdxi*invJ0;
            
            % particle mass and momentum to node
            stress = s(pid,:);
            
            for i=1:length(esctr)
                id    = esctr(i);
                dNIdx = dNdx(i,1);
                dNIdy = dNdx(i,2);
                nmass(id)       = nmass(id)       + N(i)*Mp(pid);
                nmomentum(id,:) = nmomentum(id,:) + N(i)*Mp(pid)*vp(pid,:);
                niforce(id,1)   = niforce(id,1) - Vp(pid)*(stress(1)*dNIdx + stress(3)*dNIdy);
                niforce(id,2)   = niforce(id,2) - Vp(pid)*(stress(3)*dNIdx + stress(2)*dNIdy);
            end                                                            
        end
    end
    
    % debug
            
    % update nodal momenta
    
    nmomentum = nmomentum + niforce*dtime;
 
    
    % update particle velocity and position and stresses
    k = 0;
    u = 0;
    for e=1:elemCount
        esctr = element(e,:);
        enode = node(esctr,:);
        mpts  = mpoints{e};        
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            pt(1)= (2*xp(pid,1)-(enode(1,1)+enode(2,1)))/deltax;
            pt(2)= (2*xp(pid,2)-(enode(2,2)+enode(3,2)))/deltay;
            
            [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
            J0       = enode'*dNdxi;             % element Jacobian matrix
            invJ0    = inv(J0);
            dNdx     = dNdxi*invJ0;
            Lp = zeros(2,2);
            for i=1:length(esctr)
                id = esctr(i);     
                vI = [0 0];
                if nmass(id) > tol
                    vp(pid,:)  = vp(pid,:) + dtime * N(i)*niforce(id,:)  /nmass(id);
                    xp(pid,:)  = xp(pid,:) + dtime * N(i)*nmomentum(id,:)/nmass(id);
                    vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
                end                                                
                Lp = Lp + vI'*dNdx(i,:);         % particle gradient velocity 
            end
                                  
            F       = ([1 0;0 1] + Lp*dtime)*reshape(Fp(pid,:),2,2);
            Fp(pid,:)= reshape(F,1,4);
            Vp(pid) = det(F)*Vp0(pid);                                    
            dEps    = dtime * 0.5 * (Lp+Lp');    
            dsigma  = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
            s(pid,:)  = s(pid,:) + dsigma';     
            eps(pid,:)= eps(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];     
            
            k = k + 0.5*(vp(pid,1)^2+vp(pid,2)^2)*Mp(pid);
            u = u + 0.5*Vp(pid)*s(pid,:)*eps(pid,:)';
                
%                 
%             for i=1:length(esctr)
%               id = esctr(i);                            
%               gStress(:,id)  = gStress(:,id) + N(i)*s(pid,:)';                            
%             end   
        end
    end
    
    % store time,velocty for plotting
    
    pos{istep} = xp;
    vel{istep} = vp;
    
    % update the element particle list
    
    for p=1:pCount
        x = xp(p,1);
        y = xp(p,2);
        e = floor(x/deltax) + 1 + mesh.numx*floor(y/deltay);
        pElems(p) = e;
    end
    
    for e=1:elemCount
        id  = find(pElems==e);
        mpoints{e}=id;
    end
    
    % store time,velocty for plotting
    
    ta = [ta;t];   
    ka = [ka;k];
    sa = [sa;u];
  
    % VTK output
    
    if (  mod(istep,interval) == 0 )
        xp = pos{istep};        
        data.stress  = [s zeros(pCount,1)];
        vtkFile = sprintf('../../results/mpm/twoDisks/%s%d',vtkFileName,istep);
        VTKParticles(xp,vtkFile,data);
        %vtkFileName1 = sprintf('../results/%s%d','mpmTwoDisksGrid',istep);
        %VTKPostProcess(node,element,2,'Quad4',vtkFileName1,gStress', gDisp');
    end
    
        
    % advance to the next time step
    
    t     = t + dtime;
    istep = istep + 1;
end

%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

%%


% pvdFile = fopen(strcat('../results/',vtkFileName,'.pvd'), 'wt');
% 
% fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
% fprintf(pvdFile,'<Collection>\n');
% 
% for i = 1:nsteps
%      if (  mod(i,interval) == 0 )
%     vtuFile = sprintf('%s%d%s',vtkFileName,i,'.vtp');
%     fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
%      end
% end
% 
% fprintf(pvdFile,'</Collection>\n');
% fprintf(pvdFile,'</VTKFile>\n');
% 
% fclose(pvdFile);

%ka = ka*1000;
%sa = sa*1000;

figure 
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
plot(ta(1:end),sa(1:end),'r--','LineWidth',2);
plot(ta(1:end),ka(1:end)+sa(1:end),'g-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy (x1E-3)')
legend('kinetic','strain','total')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
%axis([0 3 0 3])


% savefile = 'mpm-explicit-2disks-dt001.mat';
% save(savefile,'ta','ka','sa');

disp([num2str(toc),'   DONE '])
