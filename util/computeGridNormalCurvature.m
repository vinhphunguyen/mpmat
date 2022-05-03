function [res] = computeGridNormalCurvature(grid,body) 

% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia
% 3 September 2015.

%% retrieve grid data
node      = grid.node;
element   = grid.element;         % all grid elements
elemCount = grid.elemCount;
nodeCount = length(node);
numx2     = grid.numx;
numy2     = grid.numy;
deltax    = grid.deltax;
deltay    = grid.deltay;
% retrieve body data
elems     = body.elements;        % elements contain body under consideration
mpoints   = body.mpoints;
xp        = body.coord;
Mp        = body.mass;

option = 0;

%% cell density computation
cellDensity  = zeros(elemCount,1);
sCellDensity = zeros(elemCount,1); % smoothed cell density
volInv      = 1/(deltax*deltay);

if ( option == 1 )
for i=1:elemCount
    ie        = i;                    
    mpts      = mpoints{ie};    
    cellDensity(ie) = sum(Mp(mpts)) * volInv;
end
else
  for i=1:elemCount
    ie        = i;
    neighbors = getNeighbors(ie, numx2, numy2);
    center    = 1/4*sum( node(element(ie,:),:) );
    for in=1:length(neighbors)
      elemId = neighbors(in);
      mpts   = mpoints{elemId};
      for p=1:length(mpts)
        pid  = mpts(p);
        x    = xp(pid,:);
        [phi,~]=getQuadraticBspline2D(x-center,deltax,deltay);
        cellDensity(ie) = cellDensity(ie) + phi*Mp(pid);
      end
    end
    cellDensity(ie) = cellDensity(ie) * volInv;
  end
end



% smooth cell density

smoothVector = [1 2 1 2 4 2 1 2 1]; % for a 3x3 neighbor
sumS         = sum(smoothVector);

for i=1:elemCount
    ie              = i;
    neighbors       = getNeighbors(ie, numx2, numy2);
    dens            = cellDensity(neighbors);
    if length(dens) == 9
      sCellDensity(ie) = dot(dens,smoothVector)/sumS;
    end
end


sCellDensity = cellDensity;
%% Now, compute normals at nodes of all elements
normals = zeros(nodeCount,2);
curva   = zeros(nodeCount,1);

for iel = 1 : length(element)           % loop over all elements
    eId    = (iel);
    sctr   = element(eId,:);
    dens   = sCellDensity(eId);
    center = mean( node(sctr,:) );
    for in=1:4
        nId        = sctr(in);
        xI         = node(nId,:);
        [~,dNdx]   = getMPM2D(center-xI,deltax,deltay);
        normals(nId,:) = normals(nId,:) + dNdx*dens;
    end
end

%% now compute curvature at element centers
curvature   = zeros(elemCount,1);  % cell-centered curvature
divLength   = zeros(2,1);          % divergence of grid normal lengths

tol         = 1e-5;

for iel = 1 : length(elems)
    eId        = elems(iel);
    sctr       = element(eId,:);
    gNormals   = normals(sctr,:);       % grid normals of current element
    center     = mean( node(sctr,:) );  % center coordinate
    cNormal    = mean( gNormals );      % cell-centerd normal, non normalised
    le         = norm( cNormal );
    %le = 1;
    if ( le < tol ), le = 1; end
    cnNormal   = cNormal/le;            % cell-centerd normal, normalised    
    
    divNormal    = 0;
    divLength(:) = 0;
    for in=1:4
        ll           = norm(gNormals(in,:));
        %if ll < tol, ll = 0; end
        nId          = sctr(in);
        xI           = node(nId,:);
        [~,dNdx]     = getMPM2D(center-xI,deltax,deltay);
        
        divNormal    = divNormal + dNdx(1)*gNormals(in,1) + dNdx(2)*gNormals(in,2);
        
        divLength(1) = divLength(1) + dNdx(1)*ll;
        divLength(2) = divLength(2) + dNdx(2)*ll;
    end      
    cur = (1/le)*( cnNormal(1)*divLength(1) + cnNormal(2)*divLength(2) - divNormal) ;   
    curvature(eId) = cur;
    
    for in=1:4
        nId        = sctr(in);
        xI         = node(nId,:);
        [N,~]      = getMPM2D(center-xI,deltax,deltay);
        curva(nId,:) = curva(nId,:) + N*cur;
    end                            
end

%%

res.gNormals    = normals;
res.cDensities  = cellDensity;
res.scDensities = sCellDensity;
res.cCurvatures = curvature;
res.gCurvatures = curva;

