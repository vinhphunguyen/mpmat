% This file implements the Generalized Material Point Method.
% One dimensional problem with material points.
% The grid is one two-noded linear element.
% 
% 1D MMS problem
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% May 2014.
% Monash Uni, September 2019

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

amps   = [0.0001 0.05];
cells  = [4 5 6 7 8 9 10 12]; % choose number of elements for convenrgence study
cpGIMP  = 1;

gridSize = zeros(length(cells),1);
errMeas  = zeros(length(cells),length(amps));
totalErrMeas  = zeros(length(cells),length(amps));

for g=1:length(amps)
    G=amps(g);
    
    % function handles for MMS (manufactured solutions)
    
    mmsU = @(x,t)      G*sin(pi*x)*sin(c*pi*t);
    mmsV = @(x,t) pi*c*G*sin(pi*x)*cos(c*pi*t);
    mmsF = @(x,t) 1 + pi*G*cos(pi*x)*sin(c*pi*t);
    mmsB = @(x,t) (1/rho)*pi^2*mmsU(x,t)*( (lambda/mmsF(x,t)/mmsF(x,t))*(1-log(mmsF(x,t))) + ...
        mu*(1+1/mmsF(x,t)/mmsF(x,t)) -E );
    
    
    %%
    %  Computational grid: two-noded elements
    ppc = 2;
    for mi=1:length(cells)
        ne          = 2^cells(mi);
        gridSize(mi) = L/ne;
        
        [mesh]    = buildGrid1D(L,ne,1); % with ghost cells
        
        elemCount = mesh.elemCount;
        nodeCount = mesh.nodeCount;
        elements  = mesh.element;
        nodes     = mesh.node;
        deltax    = mesh.deltax;
        
        particles=buildParticlesGeom(mesh,ppc,rho);
        xp  = particles.xp;
        xp0  = particles.xp;
        vp  = particles.vp;
        Vp  = particles.Vp;
        Vp0 = particles.Vp0;
        Fp  = particles.Fp;
        s   = particles.s;
        eps = particles.eps;
        Mp  = particles.Mp;
        
        % initial velocities
        
        for p=1:particles.pCount
            vp(p) = mmsV(xp(p),0);
        end
        vp0=vp;
        
        %%
        hold on
        plot(nodes,zeros(nodeCount,1)+1/4,'r-s','MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',9,'LineWidth',1.1);
        plot(particles.xp,zeros(particles.pCount,1)+1/4,'bo','MarkerEdgeColor','k',...
            'MarkerFaceColor','r',...
            'MarkerSize',9,'LineWidth',1.1);
        %axis([-10 35 0 1/2])
        
        
        %%
        % data structure to store the material points for each element
        % this data structure is updated for every time step
        
        pElems  = ones(particles.pCount ,1);
        mpoints = cell (elemCount ,1);
        
        for p=1:particles.pCount
            x = xp(p);
            e = point2ElemIndex1D(x,mesh);
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
        
        %% particle size, can be updated to do cpGIMP
        
        lp = zeros(particles.pCount,1);
        lp(:) = deltax/ppc;
        lp0   = lp;
        
        %% Time loop
        
        tol = 1e-10;
        
        dtime = 0.2*deltax/c;
        time  = 0.02;
        t     = 0;
        istep = 0;
        
        nsteps = floor(time/dtime);
        err    = zeros(nsteps,1);
        ta     = 0:dtime:time;
        totalE = 0;
        
        nmass      = zeros(nodeCount,1);  % nodal mass vector
        nmomentum0 = zeros(nodeCount,1);  % nodal momentum vector
        niforce    = zeros(nodeCount,1);  % nodal internal force vector
        neforce    = zeros(nodeCount,1);  % nodal external force vector
        
        while ( t < time )
            disp(['time step ',num2str(t)])
            nmass(:)     = 0;
            nmomentum0(:)= 0;
            niforce(:)   = 0;
            neforce(:)   = 0;
            % loop over computational cells or elements
            for e=1:elemCount
                esctr = gimpElement{e};
                mpts  = mpoints{e};
                % loop over particles
                for p=1:length(mpts)
                    pid  = mpts(p);
                    lpp = lp(pid,:);
                    mm  = Mp(pid);
                    XX  = xp0(pid);
                    vol =  Vp(pid);
                    sig = s(pid);
                    velo = vp(pid);
                    % loop over nodes contribute to this particle
                    for i=1:length(esctr)
                        id    = esctr(i);
                        x     = xp(pid,:) - nodes(id,:);
                        % shape functions and first derivatives
                        [phi,dphi]    = getGIMP(x,deltax,lpp);
                        % particle mass and momentum to node
                        nmass(id)     = nmass(id)      + phi*mm;
                        nmomentum0(id)= nmomentum0(id) + phi*mm*velo;
                        niforce(id)   = niforce(id)   - vol*sig*dphi;
                        neforce(id)   = neforce(id)   + mm*phi*mmsB(XX,t);
                    end
                end
            end
            
            % update nodal momenta
            
            % be careful with left ghost cell
            
            nmomentum0(1)  = 0; nmomentum0(nodeCount)    = 0;
            nmomentum0(2)  = 0; nmomentum0(nodeCount-1)  = 0;
            
            %     niforce(nodeCount)    = 0;
            %
            %     nmomentum(nodeCount)  = 0; % Boundary conditions f1 = m1*a1, a1=0
            
            nmomentum = nmomentum0 + niforce*dtime;
            
            nmomentum(1)  = 0; nmomentum(nodeCount)    = 0;
            nmomentum(2)  = 0; nmomentum(nodeCount-1)    = 0;
            
            % update particle velocity and position and stresses
       
            for e=1:elemCount
                esctr = gimpElement{e};
                mpts  = mpoints{e};
                % loop over particles
                for p=1:length(mpts)
                    pid = mpts(p);
                    Lp  = 0;
                    xxp = xp(pid,:);
                    lpp = lp(pid,:);
                    for i=1:length(esctr)
                        id    = esctr(i);
                        x     = xxp - nodes(id,:);
                        % shape functions and first derivatives
                        [phi,dphi]=getGIMP(x,deltax,lpp);
                        
                        if nmass(id) > 0
                            vI      = nmomentum(id)/nmass(id);
                            %if ( id == 1 ), vI = - nmomentum(3)/nmass(3); end
                            vp(pid) = vp(pid) +  phi*(nmomentum(id)-nmomentum0(id))/nmass(id);
                            xp(pid) = xp(pid) + dtime * phi*vI;
                        else
                            vI = 0.;
                        end
                        % gradient velocity
                        Lp = Lp + dphi * vI;
                    end
                    F       = (1 + Lp*dtime)*Fp(pid);
                    Fp(pid) = F;
                    Vp(pid) = F*Vp0(pid);
                    if ( cpGIMP ), lp(pid) = lp0(pid) * Fp(pid); end
                    s(pid)  = lambda*log(F)/F + mu*F - mu/F; 
                end
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
            
            for p=1:particles.pCount
                xx0 = xp0(p);
                xx  = xp(p);
                up  = xx - xx0;
                uex = mmsU(xx0,t);
                dispNorm = dispNorm + Vp0(p)*(up-uex)^2;
            end
            dispNorm = sqrt(dispNorm/sum(Vp0));
            totalE   = totalE + dispNorm;
            err(istep) = dispNorm;     
        end
        %%
        disp([num2str(toc),'   DONE ']);
        errMeas(mi,g) = max(err);
        totalErrMeas(mi,g) = totalE/particles.pCount/istep;
    end
end
%%


%% convergence plot

polyfit(log(gridSize),log(errMeas(:,1)),1)
polyfit(log(gridSize),log(errMeas(:,2)),1)

loglog(gridSize,errMeas(:,1),'black*-','LineWidth',1.8)
%loglog(gridSize,errMeasMUSL(:,1),'black*-','LineWidth',1.8)
hold on
loglog(gridSize,errMeas(:,2),'reds-','LineWidth',1.8)
loglog(gridSize,totalErrMeas(:,2),'blues-','LineWidth',1.8)
%loglog(gridSize(1:6),errorLargeMUSL(1:6),'reds-','LineWidth',1.8)
xlabel('Element size')
ylabel('Error')
%legend('MPM','Exact')
set(gca,'FontSize',16)
grid on
legend('G=0.0001, e1','G=0.05, e1', 'G=0.05, RMSE')
%set(gca ,'YTickLabel',num2str(disp,1))
% axis([0 100 -0.15 0.2])



