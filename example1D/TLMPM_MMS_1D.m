% This file implements the Total Lagrangian Material Point Method (TL MPM).
% One dimensional problem with material points.
% The grid is one two-noded linear element.
% Leapfrog time integration.
%
% Convergence test of TLMPM using the Method of Manufactured Solution (MMS).
% This file contains commands to plot a convergence curve.
%
% Vinh Phu Nguyen
% Monash University, Australia
% 6 September 2019.

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
cells  = [7]; % choose number of elements for convenrgence study

gridSize = zeros(length(cells),1);
errMeas  = zeros(length(cells),length(amps));

alpha          = 1;   % blending PIC/FLIP; 1(FLIP), 0(PIC)
beta           = 1;   % 1: standard particle position update
% 0: Leroch particle position update
doubleMapping = 1;
forceAtNodes  = 1; % calculating body forces durectly at nodes of mapping from particle forces

for g=2
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
        
        [mesh]= buildGrid1D(L,ne,0);
        
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
        
        pCount = particles.pCount;
        xp   = particles.xp;  % updated particle position
        xp0  = particles.xp;  % Xp initial coordinates
        vp   = particles.vp;  % velcotites
        Vp   = particles.Vp;  % updated volume
        Vp0  = particles.Vp0; % initial volume
        Fp   = particles.Fp;  % deformation gradient matrix
        s    = particles.s;   % stress
        eps  = particles.eps; % strain (not needed here)
        Mp   = particles.Mp;  % mass, constant
        
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
        
        dN1 = -1/deltax;
        dN2 =  1/deltax;
        
        %% Time loop
        tol = 0;
        
        dtime = 0.8*deltax/c;
        time  = 0.02;
        t     = 0.;
        istep = 0;
        
        nmass     = zeros(nodeCount,1);  % nodal mass vector
        nmomentum = zeros(nodeCount,1);  % nodal momentum vector (final)
        nmomentum0= zeros(nodeCount,1);  % nodal momentum vector (begin)
        niforce   = zeros(nodeCount,1);  % nodal internal force vector
        neforce   = zeros(nodeCount,1);  % nodal external force vector
        nvelo     = zeros(nodeCount,1);
        
        % initialisation of nodal mass
        
        for e=1:elemCount
            esctr = elements(e,:);
            enode = nodes(esctr);
            mpts  = mpoints{e};
            % loop over particles
            for p=1:length(mpts)
                pid  = mpts(p);
                xx   = xp0(pid);
                mm   = Mp(pid);
                % shape functions and first derivatives
                N1  = 1 - abs(xx-enode(1))*deltaxI;
                N2  = 1 - abs(xx-enode(2))*deltaxI;
                % particle mass and momentum to node
                nmass(esctr(1))      = nmass(esctr(1))     + N1*mm;
                nmass(esctr(2))      = nmass(esctr(2))     + N2*mm;
            end
        end
        
        nsteps = floor(time/dtime);
        err    = zeros(nsteps,1);
        ta     = 0:dtime:time;
        
        while ( t < time )
            %disp(['time step ',num2str(t)]);
            % reset data
            nmomentum0(:) = 0;
            nvelo(:)      = 0;
            niforce(:)    = 0;
            neforce(:)    = 0;
            % P2G
            for e=1:elemCount
                esctr = elements(e,:);
                node1 = esctr(1);
                node2 = esctr(2);
                enode = nodes(esctr);
                mpts  = mpoints{e};
                % loop over particles
                for p=1:length(mpts)
                    pid  = mpts(p);
                    xx   = xp0(pid);
                    mm   = Mp(pid);
                    velo = vp(pid);
                    vol  = Vp0(pid);
                    stre = s(pid);
                    
                    x1   = enode(1);
                    x2   = enode(2);
                    
                    % shape functions and first derivatives
                    N1  = 1 - abs(xx-x1)*deltaxI;
                    N2  = 1 - abs(xx-x2)*deltaxI;
                    
                    nmomentum0(node1) = nmomentum0(node1) + N1*mm*velo;
                    nmomentum0(node2) = nmomentum0(node2) + N2*mm*velo;
                    
                    % internal force
                    niforce(node1) = niforce(node1) - vol*stre*dN1;
                    niforce(node2) = niforce(node2) - vol*stre*dN2;
                    
                    % external force due to manufactured body force
                    if ( forceAtNodes == 0 )
                        neforce(node1) = neforce(node1) + mm*N1*mmsB(xx,t);
                        neforce(node2) = neforce(node2) + mm*N2*mmsB(xx,t);
                    else
                        neforce(node1) =  nmass(node1) * mmsB(x1,t);
                        neforce(node2) =  nmass(node2) * mmsB(x2,t);
                    end
                end
            end
            
            % Grid update
            nforce        = niforce + neforce;
            % leapfrog integration scheme
            %if ( istep == 0 ), nforce = 0.5*nforce; end
            nmomentum     = nmomentum0 + nforce*dtime;
            % Boundary conditions
            nmomentum0(1)          = 0;
            nmomentum0(nodeCount)  = 0;
            nmomentum(1)           = 0;
            nmomentum(nodeCount)   = 0;
            
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
                        xx   = xp0(pid);
                        mm   = Mp(pid);
                        
                        N1  = 1 - abs(xx-enode(1))*deltaxI;
                        N2  = 1 - abs(xx-enode(2))*deltaxI;
                        
                        v1   = nmomentum(node1)/nmass(node1);
                        v2   = nmomentum(node2)/nmass(node2);
                        
                        % particle velocity blending PIC/FLIP
                        % alpha = 0: PIC
                        % alpha = 1: FLIP (most MPM uses this)
                        velo     = alpha     * ( vp(pid) + N1*(nmomentum(node1)-nmomentum0(node1))/nmass(node1) + ...
                                                         + N2*(nmomentum(node2)-nmomentum0(node2))/nmass(node2) ) + ...
                            (1-alpha) * ( N1*v1 + N2*v2 );
                        
                        vp(pid) = velo;
                        
                        nvelo(node1)  = nvelo(node1)  + N1*mm*velo;
                        nvelo(node2)  = nvelo(node2)  + N2*mm*velo;
                    end
                end
                nvelo            = nvelo./nmass;
                nvelo(1)         = 0;
                nvelo(nodeCount) = 0;
            else
                nvelo = nmomentum./nmass;
            end
            
            
            % G2P
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
                    xx   = xp0(pid);
                    N1   = 1 - abs(xx-enode(1))*deltaxI;
                    N2   = 1 - abs(xx-enode(2))*deltaxI;
                    
                    if ( doubleMapping == 0 )
                        vp(pid)  = alpha*vp(pid) +  N1*(v1-alpha*v10) + N2*(v2-alpha*v20);
                        % gradient velocity
                        Lp = dN1 * v1 + dN2 * v2;
                    else
                        Lp = dN1 * nvelo(node1) + dN2 * nvelo(node2);
                    end
                    
                    xp(pid)  = xp(pid) + (1-beta)*dtime*vp(pid) + beta*dtime*(N1*v1+N2*v2);
                    
                    %F       = (1 + Lp*dtime)*Fp(pid);
                    
                    F       = Fp(pid) + dtime * Lp;
                    %F       = 1 +
                    Fp(pid) = F;
                    Vp(pid) = F*Vp0(pid);
                    s(pid)  = lambda*log(F)/F + mu*F - mu/F; % Neo-Hookean
                end
            end
            
            % update the element particle list
            %   pe = floor(xp/deltax)+1;
            %   for e=1:elemCount
            %     id  = find(pe==e);
            %     mpoints{e}=id;
            %   end
            
            % advance to the next time step
            t     = t + dtime;
            istep = istep + 1;
            
            % compute error norm
            dispNorm = 0;
            Dt = 100;
            for p=1:pCount
                xx0 = xp0(p);
                xx  = xp(p);
                J   = Fp(p);
                cc  = sqrt(J*E/rho);
                Dt  = min(deltax/cc,Dt);
                up  = xx - xx0;
                uex = mmsU(xx0,t);
                dispNorm = dispNorm + Vp0(p)*(up-uex)^2;
            end
            dispNorm = sqrt(dispNorm/sum(Vp0));
            err(istep) = dispNorm;
            disp(['time step ',num2str(Dt)]);
            
        end
        %%
        disp([num2str(toc),'   DONE ']);
        errMeas(m,g) = max(err);
    end
end


%% convergence plot

polyfit(log(gridSize),log(errMeas(:,1)),1)
polyfit(log(gridSize),log(errMeas(:,2)),1)

loglog(gridSize,errMeas(:,1),'black*-','LineWidth',1.8)
%loglog(gridSize,errMeasMUSL(:,1),'black*-','LineWidth',1.8)
hold on
loglog(gridSize,errMeas(:,2),'reds-','LineWidth',1.8)
%loglog(gridSize(1:6),errorLargeMUSL(1:6),'reds-','LineWidth',1.8)
xlabel('Element size')
ylabel('Error')
%legend('MPM','Exact')
set(gca,'FontSize',16)
grid on
legend('G=0.0001','G=0.05')
%set(gca ,'YTickLabel',num2str(disp,1))
% axis([0 100 -0.15 0.2])

