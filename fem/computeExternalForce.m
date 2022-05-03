function [fext] = computeExternalForce(mesh,g,rho)

fext      = zeros(mesh.dofCount,1);

for e=1:mesh.elemCount
  esctr = mesh.element(e,:);
  enode = mesh.node(esctr,:);
  nn    = length(esctr);
  sctr(1,1:2:2*nn) = 2*esctr-1;    
  sctr(1,2:2:2*nn) = 2*esctr  ;
  % loop over Gauss point
  for p=1:length(mesh.W)
    pt       = mesh.Q(p,:);
    [N,dNdxi]= lagrange_basis(mesh.elemType,pt);   % element shape functions
    J0       = enode'*dNdxi;                  % element Jacobian matrix
    detJ     = det(J0);
    Nmat(1,1:2:2*nn)   = N;
    Nmat(2,2:2:2*nn)   = N;    
    fext(sctr) = fext(sctr) + Nmat' * [0;g] * rho * detJ * mesh.W(p);
  end
end
