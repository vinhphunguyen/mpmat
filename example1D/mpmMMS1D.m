% This file implements the Material Point Method of Sulsky 1994.
% One dimensional problem with material points.
% The grid is one two-noded linear element.
% Leapfrog time integration.
%
% Convergence test of MPM using the Method of Manufactured Solution (MMS).
% This file contains commands to plot a convergence curve.
%
% Vinh Phu Nguyen
% Monash University, Australia
% 12 September 2019.

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
E      = 1e7;               % Young modulus
nu     = 0.0;               % Poisson ratio
rho    = 1000;              % density
lambda = E*nu/(1+nu)/(1-2*nu);
mu     = E/2/(1+nu);
c      = sqrt(E/rho);
L      = 1;
amps   = [0.0001; 0.05];
cells  = [1 2 3 4 5 6 7 8 9 10]; % choose number of elements for convenrgence study

gridSize = zeros(length(cells),1);
errMeas  = zeros(length(cells),length(amps));

% blending PIC/FLIP parameter: 1(FLIP), 0(PIC)
alpha         = 1;
doubleMapping = 1;
beta          = 1; % 1: standard particle position update
% 0: Leroch particle position update

for g=1:2
    G=amps(g);
    
    % function handles for MMS (manufactured solutions)
    
    mmsU = @(x,t)      G*sin(pi*x)*sin(c*pi*t);
    mmsV = @(x,t) pi*c*G*sin(pi*x)*cos(c*pi*t);
    mmsF = @(x,t) 1 + pi*G*cos(pi*x)*sin(c*pi*t);
    mmsB = @(x,t) (1/rho)*pi^2*mmsU(x,t)*( (lambda/mmsF(x,t)/mmsF(x,t))*(1-log(mmsF(x,t))) + ...
        mu*(1+1/mmsF(x,t)/mmsF(x,t)) -E );
    
    % [X,T] = meshgrid(0:0.05:1,0:0.0006:0.02);
    % U     = G*sin(pi*X).*sin(c*pi*T);
    % surf(X,T,U)
    % xlabel('X')
    % ylabel('t')
    % zlabel('u')
    % set(gca,'FontSize',16)
    % %camlight right
    % shading interp
    % %axis equal
    
    %%
    %  Computational grid: two-noded elements
    
    for m=1:length(cells)     
        ne          = 2^cells(m);
        gridSize(m) = L/ne;
        [mesh]      = buildGrid1D(L,ne,0);
        
        elemCount = mesh.elemCount;
        nodeCount = mesh.nodeCount;
        elements  = mesh.element;
        nodes     = mesh.node;
        deltax    = mesh.deltax;
        deltaxI   = 1 / deltax;
        dN1       = -1/deltax;
        dN2       =  1/deltax;
        %%
        % Material points:
        ppc       = 2; % # particles per cell
        particles = buildParticlesGeom(mesh,ppc,rho);
        
        pCount = particles.pCount;
        xp   = particles.xp;
        xp0  = particles.xp;
        vp   = particles.vp;
        Vp   = particles.Vp;
        Vp0  = particles.Vp0;
        Fp   = particles.Fp;
        s    = particles.s;
        eps  = particles.eps;
        Mp   = particles.Mp;
        
        % initial velocities, initial stress=0
        for p=1:particles.pCount
            vp(p) = mmsV(xp(p),0);
        end
        vp0=vp;
        
        
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
        
        dtime = 0.2*deltax/c;
        time  = 0.02;
        t     = 0.;
        istep = 0;
        
        nmass     = zeros(nodeCount,1);  % nodal mass vector
        nmomentum = zeros(nodeCount,1);  % nodal momentum vector (final)
        nmomentum0= zeros(nodeCount,1);  % nodal momentum vector (begin)
        nvelo     = zeros(nodeCount,1);  %
        niforce   = zeros(nodeCount,1);  % nodal internal force vector
        neforce   = zeros(nodeCount,1);  % nodal external force vector
        
        nsteps = floor(time/dtime);
        err    = zeros(nsteps,1);
        ta     = 0:dtime:time;
        
        while ( t < time )
            disp(['time step ',num2str(t)]);
            nmass(:)     = 0;
            nmomentum0(:)= 0;
            nvelo(:)     = 0;
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
                    xx   = xp(pid);
                    mm   = Mp(pid);
                    velo = vp(pid);
                    vol  = Vp(pid);
                    sigma=s(pid);
                    % shape functions and first derivatives
                    N1  = 1 - abs(xx-enode(1))*deltaxI;
                    N2  = 1 - abs(xx-enode(2))*deltaxI;
                    % particle mass and momentum to node
                    nmass(node1)      = nmass(node1)    + N1*mm;
                    nmass(node2)      = nmass(node2)    + N2*mm;
                    nmomentum0(node1) = nmomentum0(node1) + N1*mm*velo;
                    nmomentum0(node2) = nmomentum0(node2) + N2*mm*velo;
                    % internal force
                    niforce(node1) = niforce(node1) - vol*sigma*dN1;
                    niforce(node2) = niforce(node2) - vol*sigma*dN2;
                    % external force due to manufactured body force
                    XX   = xp0(pid);
                    neforce(node1) = neforce(node1) + mm*N1*mmsB(XX,t);
                    neforce(node2) = neforce(node2) + mm*N2*mmsB(XX,t);
                end
            end
            
            % update nodal momenta
            nmomentum = nmomentum0 + niforce*dtime;
            % Dirichlet Boundary conditions
            nmomentum0(1)  = 0; nmomentum0(nodeCount)  = 0;
            nmomentum (1)  = 0; nmomentum (nodeCount)  = 0;
            
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
                nvelo(1) = 0; nvelo(nodeCount) = 0;
            else
                nvelo = nmomentum./nmass;
            end
            
            
            % update particle velocity and position and stresses (G2P)
            
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
                    F       = (1 + Lp*dtime)*Fp(pid);
                    Fp(pid) = F;
                    Vp(pid) = F*Vp0(pid);
                    s(pid)  = lambda*log(F)/F + mu*F - mu/F; % Neo-Hookean (Cauchy=1st PK stress)
                end
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
            
            % compute error norm
            dispNorm = 0;
            for p=1:pCount
                xx0 = xp0(p);
                xx  = xp(p);
                up  = xx - xx0;
                uex = mmsU(xx0,t);
                dispNorm = dispNorm + Vp(p)*(up-uex)^2;
            end
            dispNorm = sqrt(dispNorm);
            err(istep) = dispNorm;
        end
        %%
        disp([num2str(toc),'   DONE ']);
        errMeas(m,g) = max(err);
    end  % end of loop over meshes    
end      % end of loop over amplitudes

%% convergence plot

errMeasMUSL = errMeas; 
load('MMS-1D-MPM-USL.mat','errMeas')
errMeasUSL = errMeas; 

errorLargeMUSL = errMeasMUSL(:,2);
errorLargeUSL  = errMeasUSL(:,2);
% get convergence rate
polyfit(log(gridSize),log(errMeas(:,1)),1)
polyfit(log(gridSize),log(errMeas(:,2)),1)

loglog(gridSize,errMeasUSL(:,1),'black*-','LineWidth',1.8)
%loglog(gridSize,errMeasMUSL(:,1),'black*-','LineWidth',1.8)
hold on
loglog(gridSize(1:6),errorLargeUSL(1:6),'reds-','LineWidth',1.8)
%loglog(gridSize(1:6),errorLargeMUSL(1:6),'reds-','LineWidth',1.8)
xlabel('Element size')
ylabel('Error')
%legend('MPM','Exact')
set(gca,'FontSize',16)
grid on
legend('G=0.0001','G=0.05')
%set(gca ,'YTickLabel',num2str(disp,1))
% axis([0 100 -0.15 0.2])

%%
figure
plot(ta(2:end),err);
set(gca,'FontSize',16)

