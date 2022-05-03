function data = getCPDIQuadDataGeneral(pid,particle,mesh)
% Compute CPDI-Q4 data at particle "pid".
% Generalised version with particle is a cell.
% Used for sub-sampling technique.
%
% Inputs:
%
% pid:      particle index
% particle: particle mesh
% mesh:     background mesh/grid
%
% Output:
%
% data.wf:    function weights
% data.wf:    gradient weights
% data.nodes: indices of nodes influencing particle "pid"
%
% VP Nguyen
% 14 August, 2015
% The University of Adelaide

% four corners of the particle domain

nodeIds   = particle.elem{pid};
corners   = particle.node(nodeIds,:);
nodeCount = length(nodeIds);

if ( nodeCount == 4 ) % standard Q4 
  x1     = corners(1,1);  y1     = corners(1,2);
  x2     = corners(2,1);  y2     = corners(2,2);
  x3     = corners(3,1);  y3     = corners(3,2);
  x4     = corners(4,1);  y4     = corners(4,2);
  
  % particle domain area
  
  Vp     = 0.5*   ( x1*y2  - x2*y1 + x2*y3  - x3*y2 ...
                  + x3*y4  - x4*y3 + x4*y1  - x1*y4 );
  
  % function and gradient weights
  
  c1   = (x2-x1)*(y4-y1) - (y2-y1)*(x4-x1);  
  c2   = (x2-x1)*(y3-y2) - (y2-y1)*(x3-x2);  
  c3   = (x3-x4)*(y4-y1) - (y3-y4)*(x4-x1);  
  c4   = (x3-x4)*(y3-y2) - (y3-y4)*(x3-x2);

  wg    = zeros(4,2); 
  wf    = (1/(36*Vp))*[4*c1+2*c2+2*c3+c4 ...
                       2*c1+4*c2+c3+2*c4 ...
                       c1+2*c2+2*c3+4*c4 ...
                       2*c1+c2+4*c3+2*c4];
                     

                     %wf(1)
                     %(1/Vp)*[(-1/6)*(x2-x4)*y1+(1/12)*(2*x1-x3-x4)*y2+(1/12)*(x2-x4)*y3-(1/12)*(2*x1-x2-x3)*y4]
  
  wg(1,:) = [y2-y4 x4-x2];
  wg(2,:) = [y3-y1 x1-x3];
  wg(3,:) = -wg(1,:);
  wg(4,:) = -wg(2,:);
  
  wg(:,:) = (1/(2*Vp))*wg(:,:);
  
  elems  = zeros(4,1); % indices of elements of 4 corners
  
  % find elements contain the corners
  
  for c=1:4
    xc        = corners(c,:);
    elems(c)  = point2ElemIndex(xc,mesh);
  end
  
  % nodes I where phi_I(xp) are non-zero
  
  nodes = unique(mesh.element(elems,:));
  
  data.nodes = nodes;
  data.wf    = wf;
  data.wg    = wg;
  data.Vp    = Vp;
else
  nodesA = []; wfA = []; wgA = []; vpA = [];
  for i=1:2        
    nids   = nodeIds(4*i-3:4*i);
    corners= particle.node(nids,:);
    x1     = corners(1,1);  y1     = corners(1,2);
    x2     = corners(2,1);  y2     = corners(2,2);
    x3     = corners(3,1);  y3     = corners(3,2);
    x4     = corners(4,1);  y4     = corners(4,2);    
    % particle domain area    
    Vp     = 0.5*   ( x1*y2  - x2*y1 + x2*y3  - x3*y2 ...
                    + x3*y4  - x4*y3 + x4*y1  - x1*y4 );
    
    % function and gradient weights    
    c1   = (x2-x1)*(y4-y1) - (y2-y1)*(x4-x1);
    c2   = (x2-x1)*(y3-y2) - (y2-y1)*(x3-x2);
    c3   = (x3-x4)*(y4-y1) - (y3-y4)*(x4-x1);
    c4   = (x3-x4)*(y3-y2) - (y3-y4)*(x3-x2);
    
    wg    = zeros(4,2);
    wf    = (1/(36*Vp))* [4*c1+2*c2+2*c3+c4 2*c1+4*c2+c3+2*c4 c1+2*c2+2*c3+4*c4 2*c1+c2+4*c3+2*c4];
    
    wg(1,:) = [y2-y4 x4-x2];
    wg(2,:) = [y3-y1 x1-x3];
    wg(3,:) = -wg(1,:);
    wg(4,:) = -wg(2,:);
    
    wg(:,:) = (1/(2*Vp))*wg(:,:);
    
    elems  = zeros(4,1); % indices of elements of 4 corners
    
    % find elements contain the corners
    
    for c=1:4
      xc        = corners(c,:);
      elems(c)  = point2ElemIndex(xc,mesh);
    end
    
    % nodes I where phi_I(xp) are non-zero
    
    nodes = mesh.element(elems,:);

    nodesA = [nodesA;nodes];
    wfA    = [wfA;wf'];
    wgA    = [wgA;wg];
    vpA    = [vpA;Vp];
  end
      
  data.nodes = unique(nodesA);
  data.wf    = wfA;
  data.wg    = wgA;
  data.Vp    = vpA;
end


