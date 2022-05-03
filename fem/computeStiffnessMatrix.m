function [stiffMat] = computeStiffnessMatrix(mesh,material)
%
% This function computes the stiffness matrix for a FE mesh stored in 'mesh'. 
% Only 2D problems.
% Vinh Phu Nguyen

stiffMat  = zeros(mesh.dofCount,mesh.dofCount);

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
    dNdx     = dNdxi/J0;  
    B(1,1:2:2*nn)  = dNdx(:,1)';
    B(2,2:2:2*nn)  = dNdx(:,2)';
    B(3,1:2:2*nn)  = dNdx(:,2)';
    B(3,2:2:2*nn)  = dNdx(:,1)';    
    stiffMat(sctr,sctr)  = stiffMat(sctr,sctr) + B'*material.D*B*detJ*mesh.W(p);    
  end
end
