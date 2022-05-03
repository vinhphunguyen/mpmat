function [particles] = updateCorners (pCorners,particles,mesh,nmomentum,nmass,dtime)

for c=1:length(pCorners)              % loop over corners of 'p'
  cid   = pCorners(c);                % index of corner 'c'
  xc    = particles.node(cid,:);      % coordinates of corner 'c'
  ec    = point2ElemIndex(xc,mesh);   % element contains corner 'c'
  esctr = mesh.element(ec,:);
  for i=1:length(esctr)
    nid       = esctr(i);
    x         = xc - mesh.node(nid,:);
    [N,~]     = getMPM2D(x,mesh.deltax,mesh.deltay);
    if nmass(nid) > 0
      xc = xc + dtime*N*nmomentum(nid,:)/nmass(nid);
    end
  end
  particles.node(cid,:) = xc;
end


