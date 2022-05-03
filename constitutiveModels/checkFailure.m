function [isFailure] = checkFailure(ai)

global crackedElems notchedElems cracks mat a
global mesh quad u0 stress0 stress 
global status normal

disp('checking failure...')

isFailure=0;

for e=1:mesh.elemCount    
  %if (ismember(e,notchedElems)) continue;end
  
  if ~isempty(notchedElems)
    if (ismember(e,notchedElems)) continue;end
  end
  
  % the following is to allow cracking for only some elements
  if ~isempty(crackedElems)
    if (~ismember(e,crackedElems)) continue;end
  end
  
  for q=1:size(quad.W,1)                          % quadrature loop
    if status(q,e)==1 continue; end  
    sigma  = stress(:,q,e);
    sigmai = getPrincipalStress(sigma,mat.nu,mat.stressState);
    yield  = sigmai(1) - mat.ft0;
    ft     = sigmai(1);
    if (yield>=0)      
      isFailure     = 1;
      status(q,e)   = 1;
      mat.fts(q,e)  = ft;
      ai(:,q,e)     = [-ft^2/mat.Gf 0 0 mat.ks];
      n             = getPrincipalDirection( sigma );
      %n=[1 0];
      normal(:,q,e) = n;
      % compute a crack segment, for visualisation purpose
      pt       = quad.Q(q,:);                       % quadrature point
      wt       = quad.W(q);                         % quadrature weight
      [N,~]= lagrange_basis(mesh.elemType,pt);
      sctr = mesh.element(e,:);
      xp   = N'*mesh.node(sctr,:);
      tem1 = 0.25*mesh.H(e)*n(2);
      tem2 = 0.25*mesh.H(e)*n(1);
      cracks(:,q,e) = [xp(1)-tem1 xp(2)-tem2 xp(1)+tem1 xp(2)+tem2];
    end
  end                                             % of quadrature loop
end                                               % of element loop





