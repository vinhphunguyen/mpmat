function particles=buildParticles(mesh,ppc,rho)

% mesh: Eulerian grid
% ppc : number of particles per cell

if (ppc > 1)
    [W,Q]=quadrature( ppc, 'GAUSS', 1 );
else
    W = 2;
    Q = 0;
end

Mp   = [];
Vp   = [];
xp   = [];

for e=1:mesh.elemCount                 % start of element loop
    if (mesh.ghostCell)
        if (e==1) || (e==mesh.elemCount)
            continue;
        end
    end
    sctr = mesh.element(e,:);          %  element scatter vector
    pts  = mesh.node(sctr,:);
    for q=1:size(W,1)                           % quadrature loop
        pt=Q(q,:);                              % quadrature point
        wt=W(q);                                % quadrature weight
        [N,dNdxi]=lagrange_basis('L2',pt);
        J0 = norm(pts'*dNdxi);
        x  = N'*pts;
        
        Vp  = [Vp;wt*J0];
        Mp  = [Mp;wt*J0*rho];
        xp  = [xp;x];
    end
end

pCount = length(xp);

particles.xp  = xp;
particles.Mp  = Mp;
particles.Vp  = Vp;                            % volume
particles.Vp0 = Vp;
particles.Fp  = ones(pCount,1);                % gradient deformation
particles.s   = zeros(pCount,1);               % stress
particles.eps = zeros(pCount,1);               % strain
particles.vp  = zeros(pCount,1);               % velocity
particles.pCount=pCount;
