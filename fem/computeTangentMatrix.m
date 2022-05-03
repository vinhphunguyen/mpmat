function [geoMat, stiffMat, fint] = computeTangentMatrix(mesh,material,ndisp)
%
% This function computes the internal force vector, material tangent matrix
% and geometry tangent matrix for a FE mesh stored in 'mesh'. 
% Only for neo-Hookean materials.
% Only 2D problems.
% Vinh Phu Nguyen

identity  = [1 0;0 1]; 

fint      = zeros(mesh.dofCount,1);
stiffMat  = zeros(mesh.dofCount,mesh.dofCount);
geoMat    = zeros(mesh.dofCount,mesh.dofCount);
H         = zeros(mesh.nnode,mesh.nnode);
Kgeo      = zeros(mesh.ndof ,mesh.ndof);

for e=1:mesh.elemCount
  esctr = mesh.element(e,:);
  enode = mesh.node(esctr,:);
  nn    = length(esctr);
  sctr(1,1:2:2*nn) = 2*esctr-1;    
  sctr(1,2:2:2*nn) = 2*esctr  ;
  ue               = [ndisp(2*esctr-1)'; ndisp(2*esctr)'];
  H(:) = 0;
  % loop over Gauss point
  for p=1:length(mesh.W)
    pt       = mesh.Q(p,:);
    [N,dNdxi]= lagrange_basis(mesh.elemType,pt);   % element shape functions
    J0       = enode'*dNdxi;                  % element Jacobian matrix
    detJ     = det(J0);
    dNdx     = dNdxi/J0;  
    % compute stress (Cauchy) at this GP
    F        = inv(identity - ue*dNdx);
    %F        = enode'*dNdx;
    invF     = inv(F);
    detF     = det(F);
    %if detF < 0, error('negative J'); end
    sig      = (material.mu*(F*F'-identity) + material.lambda*log(detF)*identity)/detF;
    % material tangent matrix D 
    mu       = material.mu - material.lambda*log(detF);
    D(1,1)   = material.lambda + 2*mu; D(1,2)   = material.lambda;
    D(2,1)   = material.lambda; D(2,2) = D(1,1);
    D(3,3)   = mu;
    B(1,1:2:2*nn)  = dNdx(:,1)';
    B(2,2:2:2*nn)  = dNdx(:,2)';
    B(3,1:2:2*nn)  = dNdx(:,2)';
    B(3,2:2:2*nn)  = dNdx(:,1)';
    fint(sctr)           = fint(sctr)          + B'*[sig(1,1);sig(2,2);sig(1,2)]*detJ*mesh.W(p);
    stiffMat(sctr,sctr)  = stiffMat(sctr,sctr) + B'*D*B*detJ*mesh.W(p);
    H                    = H + dNdx*sig*dNdx'*detJ*mesh.W(p);
  end
  
  Kgeo(1,1:2:2*nn) = H(1,:);
  Kgeo(1,2:2:2*nn) = H(1,:);
  Kgeo(2,1:2:2*nn) = H(2,:);
  Kgeo(2,2:2:2*nn) = H(2,:);
  Kgeo(3,1:2:2*nn) = H(3,:);
  Kgeo(3,2:2:2*nn) = H(3,:);
  Kgeo(4,1:2:2*nn) = H(4,:);
  Kgeo(4,2:2:2*nn) = H(4,:);
  
  geoMat(sctr,sctr)  = geoMat(sctr,sctr) + Kgeo;
end
