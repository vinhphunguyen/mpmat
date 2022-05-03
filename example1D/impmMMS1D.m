% This file implements the Improved Material Point Method (Std MPM with MLS data reconstruction).
% One dimensional problem with material points.
% The grid is one two-noded linear element.
% Leapfrog time integration.
%
% Moving Least Square approximation used to construct particle data on the
% grid (nodes/centers).
% MLS introduces two parameters into the problem: smoothing length (domain of influence).
% Generally, vvelocity MLS has one smoothing length and density/stress MLS have
% another smoothing length. Parameter study is needed.
%
% Convergence test of MPM using the Method of Manufactured Solution (MMS).
% This file contains commands to plot a convergence curve.
%
% Vinh Phu Nguyen
% Monash University, Australia
% 16 Sep 2019.

%%

addpath ../grid/
addpath ../basis/
addpath ../particleGen/
addpath ../constitutiveModels/
addpath ../util/
addpath ../geoMesh/
addpath ../externals/
addpath ../postProcessing/
addpath ../mls/

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
cells  = [4 5 6 7 8 9 10 12]; % choose number of elements for convenrgence study

gridSize = zeros(length(cells),1);
errMeas  = zeros(length(cells),length(amps));

forceAtGrid = 0;

for g=1:length(amps)
    G=amps(g);
    
    % function handles for MMS (manufactured solutions)
    
    mmsU = @(x,t)      G*sin(pi*x)*sin(c*pi*t);
    mmsV = @(x,t) pi*c*G*sin(pi*x)*cos(c*pi*t);
    mmsF = @(x,t) 1 + pi*G*cos(pi*x)*sin(c*pi*t);
    mmsB = @(x,t) (1/rho)*pi^2*mmsU(x,t)*( (lambda/mmsF(x,t)/mmsF(x,t))*(1-log(mmsF(x,t))) + ...
        mu*(1+1/mmsF(x,t)/mmsF(x,t)) -E );
    
    %% MLS weight functions
    shape = 'circle' ;         % shape of domain of influence
    dmax1 = 2.0 ;              % radius = dmax * nodal spacing
    dmax2 = 2.0 ;              % radius = dmax * nodal spacing
    form  = 'cubic_spline' ;   % using cubic spline weight function
    
    %%
    %  Computational grid: two-noded elements
    
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
        omegaC    = deltax;
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
        rhop = particles.rho;
        
        % initial velocities, initial stress=0
        for p=1:particles.pCount
            vp(p) = mmsV(xp(p),0);
        end
        vp0=vp;
        
        % Domain of influence for every nodes(sampling points)
        % Uniformly distributed nodes
        % Definition : rad = dmax*deltaX
        
        di1 = ones(pCount,1)*dmax1*deltax;
        di2 = ones(pCount,1)*dmax2*deltax;
        
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
        
        %% nodal quantities
        nmass       = zeros(nodeCount,1);  % nodal mass vector
        nvelo       = zeros(nodeCount,1);  % nodal velocity vector (final)
        nvelo0      = zeros(nodeCount,1);  % nodal velocity vector (begin)
        niforce     = zeros(nodeCount,1);  % nodal internal force vector
        neforce     = zeros(nodeCount,1);  % nodal external force vector
        
        cellDensity = zeros(elemCount,1);  % cell-centerd density
        cellStress  = zeros(elemCount,1);  % cell-centerd stress
        
        for e=1:elemCount
            esctr = elements(e,:);
            enode = nodes(esctr);
            xc    = mean(enode);
            index = defineSupport(xp,xc,di2);
            %if length(index) <= 2, disp('A singular'); end
            phi   = mlsLinearBasis1D(xc,index,xp,di2,form);
            cellDensity(e) = cellDensity(e) + dot(phi,rhop(index));
            cellStress(e)  = cellStress(e)  + dot(phi,s(index));
        end
        
        
        %% Time loop
        tol = 0;
        
        dtime = 0.2*deltax/c;
        time  = 0.02;
        t     = 0.;
        istep = 0;
        
        nsteps = floor(time/dtime);
        err    = zeros(nsteps,1);
        ta     = 0:dtime:time;
        
        while ( t < time )
            disp(['time step ',num2str(t)]);
            nmass(:)     = 0;
            nvelo0(:)    = 0;
            niforce(:)   = 0;
            neforce(:)   = 0;
            % loop over computational cells or elements
            for e=1:elemCount
                esctr = elements(e,:);
                node1 = esctr(1);
                node2 = esctr(2);
                enode = nodes(esctr);
                mpts  = mpoints{e};
                % one point quadrature
                xc   = mean(enode);
                % shape functions and first derivatives
                N1  = 0.5; %1 - abs(xc-enode(1))*deltaxI;
                N2  = 0.5; %1 - abs(xc-enode(2))*deltaxI;
               
                % mass
                rhoC           = cellDensity(e);
                nmass(node1)   = nmass(node1)     + N1*rhoC*omegaC;
                nmass(node2)   = nmass(node2)     + N2*rhoC*omegaC;
                % internal force
                sigma          = cellStress(e);
                niforce(node1) = niforce(node1) - omegaC*sigma*dN1;
                niforce(node2) = niforce(node2) - omegaC*sigma*dN2;
                % external force due to manufactured body force
                % there are 3 ways:
                % 1. from cell centered similar to mass and internal force
                % 2. from particles as in standard MPM (better than 1)
                % 3. directly at the nodes
             
                if ( forceAtGrid == 0 )
                    for p=1:length(mpts)
                        pid  = mpts(p);
                        xx   = xp(pid);
                        mm   = Mp(pid);
                        % shape functions and first derivatives
                        N1  = 1 - abs(xx-enode(1))*deltaxI;
                        N2  = 1 - abs(xx-enode(2))*deltaxI;
                        % external force due to manufactured body force
                        XX   = xp0(pid);
                        neforce(node1) = neforce(node1) + mm*N1*mmsB(XX,t);
                        neforce(node2) = neforce(node2) + mm*N2*mmsB(XX,t);
                    end
                else
                      neforce(node1) = neforce(node1) + rhoC*omegaC*N1*mmsB(xc,t);
                      neforce(node2) = neforce(node2) + rhoC*omegaC*N2*mmsB(xc,t);
%                      neforce(node1) =  nmass(node1) * mmsB(x1,t);
%                      neforce(node2) =  nmass(node2) * mmsB(x2,t);
                end
            end
            % project velocity to grid nodes (MLS)
            for i=1:nodeCount
                pt    = nodes(i);
                index = defineSupport(xp,pt,di1);
                %if length(index) <= 2, disp('A singular'); end
                phi   = mlsLinearBasis1D(pt,index,xp,di1,form);
                nvelo0(i) = nvelo0(i) + dot(phi,vp(index));
            end
            
            % update nodal velocity
            nforce        = niforce + neforce;
            % leapfrog integration scheme
            %if ( istep == 0 ), nforce = 0.5*nforce; end
            nvelo     = nvelo0 + dtime*nforce./nmass;
            % Boundary conditions
            nvelo(1)  = 0;nvelo0(1)  = 0;
            nvelo(nodeCount)  = 0; nvelo0(nodeCount)  = 0; 
            % update particle velocity and position and stresses
            for e=1:elemCount
                esctr = elements(e,:);
                node1 = esctr(1);
                node2 = esctr(2);
                enode = nodes(esctr);
                mpts  = mpoints{e};
                % loop over particles
                for p=1:length(mpts)
                    pid = mpts(p);
                    xx  = xp(pid);
                    N1  = 1 - abs(xx-enode(1))*deltaxI;
                    N2  = 1 - abs(xx-enode(2))*deltaxI;
                    dN1 = -1/deltax;
                    dN2 =  1/deltax;
                    
                    v1       = nvelo(node1);
                    v2       = nvelo(node2);
                   
                    vp(pid)  = vp(pid) + N1*(v1-nvelo0(node1)) + ...
                                       + N2*(v2-nvelo0(node2));

                    xp(pid)  = xp(pid) + dtime * ( N1*v1 +  N2*v2 );
                    

                    % gradient velocity
                    Lp      = dN1 * v1 + dN2 * v2;
                    F       = (1 + Lp*dtime)*Fp(pid);
                    Fp(pid) = F;
                    Vp(pid) = F*Vp0(pid);
                    s(pid)  = lambda*log(F)/F + mu*F - mu/F; % Neo-Hookean
                    rhop(pid) = rho/F;
                end
            end
            
            % project particle density/stress to grid centers (MLS)
            cellDensity(:) = 0;
            cellStress(:)  = 0;
            for e=1:elemCount
                esctr = elements(e,:);
                enode = nodes(esctr);
                xc    = mean(enode);
                index = defineSupport(xp,xc,di2);
                %if length(index) <= 2, disp('A singular'); end
                phi   = mlsLinearBasis1D(xc,index,xp,di2,form);
                cellDensity(e) = cellDensity(e) + dot(phi,rhop(index));
                cellStress(e)  = cellStress(e)  + dot(phi,s(index));
            end
            
            % update the element particle list
            pe = floor(xp/deltax)+1;
            for e=1:elemCount
                id  = find(pe==e);
                mpoints{e}=id;
            end
            
            % advance to the next time step
            t     = t + dtime;
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
        end % end of time loop 
        errMeas(m,g) = max(err);
    end     % end of mesh loop
end         % end of G loop
%%
disp([num2str(toc),'   DONE ']);


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


%%
figure
plot(ta(2:end),err);
set(gca,'FontSize',16)


