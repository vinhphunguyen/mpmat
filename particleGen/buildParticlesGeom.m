function particles=buildParticlesGeom(mesh,ppc,rho)
% distribute particles in a 1D structured mesh. 
% Using geometry not determined by Gauss point position
% as the function buildParticles(...).
% Inputs:
% mesh: Eulerian grid
% ppc : number of particles per cell
% 
% VP Nguyen
% May 2014, Saigon, Vietnam

Mp   = [];
Vp   = [];
xp   = [];

dx   = mesh.deltax/(ppc+1);
vol  = mesh.deltax/(ppc);

for e=1:mesh.elemCount                 % start of element loop
    if (mesh.ghostCell)
        if (e==1) || (e==mesh.elemCount)
            continue;
        end
    end
    sctr = mesh.element(e,:);          %  element scatter vector
    pts  = mesh.node(sctr,:);
    for q=1:ppc                                
        %x   = pts(1) + q*dx;        
        x   = pts(1) + vol*0.5 + (q-1)*vol;        
        Vp  = [Vp;vol];
        Mp  = [Mp;vol*rho];
        xp  = [xp; x];
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
particles.T   = zeros(pCount,1);               % temperature
particles.C   = zeros(pCount,1);               % specific heat
particles.q   = zeros(pCount,1);               % heat flux q
particles.pCount=pCount;
particles.rho = rho*ones(pCount,1);
