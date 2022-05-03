% This file implements the Material Point Method of Sulsky 1994.
% With 1D CPDI of Brannon.
% One dimensional problem with material points.
% The grid is one two-noded linear element.
%
% Convergence test of CPDI using the Method of Manufactured Solution (MMS).
% This file contains commands to plot a convergence curve.
% STATUS:
% Cannot obtain convergence of CPDI for this (MPM can!!!). For the
% vibration bar, CPDI does converge at 1st order!!! THere something wrong
% in either my implementation or CPDI theory.
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
E      = 1e7;               % Young modulus
nu     = 0.0;               % Poisson ratio
rho    = 1000;              % density
lambda = E*nu/(1+nu)/(1-2*nu);
mu     = E/2/(1+nu);
c      = sqrt(E/rho);
L      = 1;
amps   = [0.0001; 0.05];
cells  = [4]; % choose number of elements for convenrgence study

gridSize = zeros(length(cells),1);
errMeas  = zeros(length(cells),length(amps));

for g=1:2
    G=amps(g);
    
    % function handles for MMS (manufactured solutions)
    
    mmsU = @(x,t)      G*sin(pi*x)*sin(c*pi*t);
    mmsV = @(x,t) pi*c*G*sin(pi*x)*cos(c*pi*t);
    mmsF = @(x,t) 1 + pi*G*cos(pi*x)*sin(c*pi*t);
    mmsB = @(x,t) (1/rho)*pi^2*mmsU(x,t)*( (lambda/mmsF(x,t)/mmsF(x,t))*(1-log(mmsF(x,t))) + ...
        mu*(1+1/mmsF(x,t)/mmsF(x,t)) -E );
    
    
    %%
    %  Computational grid: two-noded elements
    
    for mi=1:length(cells)
        ne          = 2^cells(mi);
        gridSize(mi) = L/ne;
        np          = ne*2;
        
        [mesh]     = buildGrid1D(L,ne,1);
        
        nodeCount = mesh.nodeCount;
        elements  = mesh.element;
        nodes     = mesh.node;
        deltax    = mesh.deltax;
        
        %%
        % Material points:
        
        [particles] = buildGrid1D(L,np,0);
        if mesh.ghostCell, particles.node = particles.node + mesh.deltax; end
        particles.pCount = np;
        xp               = zeros(particles.pCount,1);
        vp               = zeros(particles.pCount,1);
        Mp               = zeros(particles.pCount,1);
        Vp               = zeros(particles.pCount,1);
        Fp               = ones (particles.pCount,1);  % be careful with F=1!!!
        s                = zeros(particles.pCount,1);  % stress
        eps              = zeros(particles.pCount,1);  % strain
        
        % particle data (mass, ...)
        
        for p=1:particles.pCount
            xp(p)   =  mean(particles.node(particles.element(p,:)));
            vp(p)   =  mmsV(xp(p),0);
            corners = particles.node(particles.element(p,:));
            Vp(p)   = corners(2)-corners(1);
            Mp(p)   = rho*Vp(p);
        end
        vp0 = vp;
        Vp0 = Vp;
        xp0 = xp;
        %q  = Mp*vp;                          % momentum
        
        %%
%         figure(1)
%         hold on
%         plot(nodes,zeros(nodeCount,1)+1/2,'r-s','MarkerEdgeColor','k',...
%             'MarkerFaceColor','g',...
%             'MarkerSize',9,'LineWidth',1.1);
%         plot(xp,zeros(particles.pCount,1)+1/2,'b*');
%         plot(particles.node,zeros(particles.pCount+1,1)+1/2,'c*','MarkerSize',16);
%         box on
%         %axis([0 L+0.5 0 1.])
        
        
        %% Time loop
        
        tol = 1e-12;
        
        dtime = 0.2*deltax/c;
        time  = 0.02;
        t     = 0;
        
        istep = 0;
        nsteps = floor(time/dtime);
        err    = zeros(nsteps,1);
        ta     = 0:dtime:time;
        
        nmass      = zeros(nodeCount,1);  % nodal mass vector
        nmomentum  = zeros(nodeCount,1);  % nodal momentum vector
        nmomentum0 = zeros(nodeCount,1);  % nodal momentum vector
        niforce    = zeros(nodeCount,1);  % nodal internal force vector
        neforce    = zeros(nodeCount,1);  % nodal external force vector
        
        while ( t < time )
            nmass(:)     = 0;
            nmomentum0(:) = 0;
            niforce(:)   = 0;
            neforce(:)   = 0;
            % loop over particles (P2G)
            for p=1:particles.pCount
                pData  = getCPDI1D(p,particles,mesh);
                inodes = pData.node;
                XX     = xp0(p);
                mp     = Mp(p);
                velo   = vp(p);
                vol    = Vp(p);
                stre   = s(p);
                for i=1:length(inodes)
                    nI = inodes(i);
                    Ni = pData.phi(i);
                    
                    % particle mass and momentum to node
                    nmass(nI)      = nmass(nI)      + Ni*mp;
                    nmomentum0(nI) = nmomentum0(nI) + Ni*mp*velo;
                    % internal force
                    niforce(nI) = niforce(nI) - vol*stre*pData.dphi(i);
                    % external force due to manufactured body force
                    neforce(nI) = neforce(nI) + mp*Ni*mmsB(XX,t);
                end
            end
            
            % update nodal momenta
            nforce        = niforce + neforce;
            % leapfrog integration scheme
            %if ( istep == 0 ), nforce = 0.5*nforce; end
            nmomentum     = nmomentum0 + nforce*dtime;
            
             % Boundary conditions f1 = m1*a1, a1=0
            nmomentum0([1,2])  = 0; nmomentum0([nodeCount,nodeCount-1])  = 0;
            nmomentum([1,2])   = 0; nmomentum([nodeCount,nodeCount-1])   = 0;
            
            % update particle velocity and stresses (G2P)
            for p=1:particles.pCount
                pData  = getCPDI1D(p,particles,mesh);
                inodes = pData.node;
                Lp     = 0;
                for i=1:length(inodes)
                    nI = inodes(i); vI = 0;
                    if nmass(nI) > tol
                        vp(p)  = vp(p) + pData.phi(i)*(nmomentum(nI)-nmomentum0(nI))/nmass(nI);
                        vI = nmomentum(nI)/nmass(nI);
                    end
                    Lp = Lp + pData.dphi(i) * vI;
                end
                F       = (1 + Lp*dtime)*Fp(p);
                Fp(p)   = F;
                Vp(p)   = F*Vp0(p);
                s(p)    = lambda*log(F)/F + mu*F - mu/F; % Neo-Hookean
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
                        %xp(pid)  = xp(pid) + dtime*(N1*v1+N2*v2);
                    end
                end
                particles.node(c) = xc+uc;
            end
            
            % advance to the next time step
            t     = t + dtime;
            istep = istep + 1;
            
            
            % compute error norm
            dispNorm = 0;
            for p=1:particles.pCount
                xx0 = xp0(p);
                xp(p) = mean(particles.node(particles.element(p,:)));
                xx  = xp(p);
                up  = xx - xx0;
                uex = mmsU(xx0,t);
                dispNorm = dispNorm + Vp(p)*(up-uex)^2;
            end
            dispNorm = sqrt(dispNorm);
            err(istep) = dispNorm;
        end % end of time loop
        errMeas(mi,g) = max(err);
        
    end     % end of meshes loop
end         % end of g loop
%%

%% convergence plot

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

