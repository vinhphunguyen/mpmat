% This file implements the Material Point Method of Sulsky 1994.
% With 1D CPDI of Brannon.
% One dimensional problem with material points.
% The grid is one two-noded linear element.
% This example is taken from "Caveats on the Implementation of the Generalized
% Material Point Method", Buzzi et al, CMES 2008.
%
% New: add convergence plot
%
% Vinh Phu Nguyen
% Monash University
% July 2016.

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

%%
%  Computational grid: two-noded elements

m             = [3 4 5 6 7];   % choose number of elements for convenrgence study
elemCount     = 2^m(3); % number of elements of the background grid
particleCount = elemCount * 1;% particle 'elements', 2 = ppc
%particleCount = 50;

[mesh]    =buildGrid1D(L,elemCount,1); % 1 => ghost cell

nodeCount = mesh.nodeCount;
elements  = mesh.element;
nodes     = mesh.node;
deltax    = mesh.deltax;

%%
% Material points:

[particles] = buildGrid1D(L,particleCount,0);

if mesh.ghostCell, particles.node = particles.node + mesh.deltax; end
particles.pCount = particleCount;
xp               = zeros(particles.pCount,1);
vp               = zeros(particles.pCount,1);
Mp               = zeros(particles.pCount,1);
Vp               = zeros(particles.pCount,1);
Fp               = ones (particles.pCount,1);  % be careful with F=1!!!
s                = zeros(particles.pCount,1);  % stress
eps              = zeros(particles.pCount,1);  % strain

% particle data (mass, ...)

for p=1:particles.pCount
    xp(p)   = mean(particles.node(particles.element(p,:)));
    vp(p)   = v0*sin(beta1*xp(p));
    corners = particles.node(particles.element(p,:));
    Vp(p)   = corners(2)-corners(1);
    Mp(p)   = rho*Vp(p);
end
vp0 = vp;
Vp0 = Vp;
xp0 = xp;
%q  = Mp*vp;                          % momentum

%%
figure(1)
hold on
plot(nodes,zeros(nodeCount,1)+1/2,'r-s','MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',9,'LineWidth',1.1);
plot(xp,zeros(particles.pCount,1)+1/2,'b*');
plot(particles.node,zeros(particles.pCount+1,1)+1/2,'c*','MarkerSize',16);
box on
%axis([0 L+0.5 0 1.])


%% Time loop

tol = 1e-5;

dtime = 0.2*deltax/c;
time  = (2*pi/omega1)*10; %time=10*dtime;
t     = 0;

ta = [];           % time
va = [];           % velocities
xa = [];           % position
ka = [];           % kinetic energy
sa = [];           % strain energy
er = [];           % errors

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,1);  % nodal momentum vector
niforce   = zeros(nodeCount,1);  % nodal internal force vector
neforce   = zeros(nodeCount,1);  % nodal external force vector

while ( t < time )
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    % loop over particles
    for p=1:particles.pCount
        pData  = getCPDI1D(p,particles,mesh);
        inodes = pData.node;
        for i=1:length(inodes)
            nI = inodes(i);
            % particle mass and momentum to node
            nmass(nI)      = nmass(nI)      + pData.phi(i)*Mp(p);
            nmomentum(nI)  = nmomentum(nI)  + pData.phi(i)*Mp(p)*vp(p);
            % internal force
            niforce(nI) = niforce(nI) - Vp(p)*s(p)*pData.dphi(i);
        end
    end
    
    % update nodal momenta
    
    nmomentum(mesh.lNodes+1)  = 0; % Boundary conditions f1 = m1*a1, a1=0
    niforce(mesh.lNodes+1)    = 0;
     %if ( t == 0 ), niforce = 0.5*niforce; end
    nmomentum = nmomentum + niforce*dtime;
    
    % update particle velocity and stresses
    k = 0;
    u = 0;
    % loop over particles
    for p=1:particles.pCount
        pData  = getCPDI1D(p,particles,mesh);
        inodes = pData.node;
        Lp     = 0;
        for i=1:length(inodes)
            nI = inodes(i); 
            vI = 0;
            if nmass(nI) > tol
                vp(p)  = vp(p) + dtime * pData.phi(i)*niforce(nI)/nmass(nI);
                vI     = nmomentum(nI)/nmass(nI);
            end            
            Lp = Lp + pData.dphi(i) * vI;
        end
        Fp(p)   = (1 + Lp*dtime)*Fp(p);
        %Vp(p)   = Fp(p)*Vp0(p);
        dEps    = dtime * Lp;
        s(p)    = s(p)   + E * dEps;
        eps(p)  = eps(p) + dEps;
        
        k = k + 0.5*vp(p)^2*Mp(p);
        u = u + 0.5*s(p)*eps(p)*Vp(p);
    end
    
    % update particle corners        
    for c=1:size(particles.node)
        xc    = particles.node(c);
        ec    = point2ElemIndex1D(xc,mesh);
        esctr = mesh.element(ec,:);
        uc    = 0;
        for i=1:length(esctr)
            id      = esctr(i);
            x       = xc - mesh.node(id,:);
            [N,dNdx]= getMPM(x,mesh.deltax);
            if nmass(id) > tol
                uc = uc + dtime*N*nmomentum(id,:)/nmass(id);            
            end
        end
        particles.node(c) = xc + uc;
    end
    
    if isnan(particles.node)
        nmass
        error('NAN in particle corners'); 
    end
% no difference compared with conventional F_p update
    for p=1:particles.pCount  
    corners = particles.node(particles.element(p,:));
    Vp(p)   = corners(2)-corners(1);
    end

                
    % advance to the next time step
    t = t + dtime;
    % store time,velocty for plotting
    cv  = 1/sum(Mp)*(dot(Mp,vp));
    
    ta = [ta;t];
    va = [va;cv];
    ka = [ka;k];
    sa = [sa;u];
    
    dispNorm = 0;
    for p=1:particles.pCount
        xx0 = xp0(p);
        xp(p) = mean(particles.node(particles.element(p,:)));
        xx  = xp(p);
        up  = xx - xx0;
        uex = (v0/omega1)*sin(omega1*t)*sin(beta1*xp(p));
        dispNorm = dispNorm + Vp(p)*(up-uex)^2;
    end
    dispNorm = sqrt(dispNorm);
    er = [er;dispNorm];
    

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
%% center of mass velocity
% disp=[1.121301946778948e-03;
%     5.937287685533104e-04;
%     3.047997447310309e-04;
%     1.543399119279315e-04];
%     %0.002400440903296];
    % velocity for all points
    disp=[1.036172019335976e-02;
    5.253431836748625e-03;
    2.640344056798967e-03;
    1.543399119279315e-04];
    %0.002400440903296];

size=[0.125000000000000;
    0.062500000000000;
    0.031250000000000;
    0.015625000000000];
    %0.007812500000000];

loglog(size,disp,'black*-','LineWidth',1.8)
hold on
xlabel('Element size')
ylabel('Error')
%legend('MPM','Exact')
set(gca,'FontSize',16)
grid on
legend('G=0.0001','G=0.01')

polyfit(log(size),log(disp),1)

%%
figure
plot(ta(1:end),er);
set(gca,'FontSize',16)

