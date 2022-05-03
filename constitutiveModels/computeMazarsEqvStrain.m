function eqvStrain = computeMazarsEqvStrain(strainl,mat)

% form strain tensor
strtens=zeros(3,3);
strtens(1,1)=strainl(1);
strtens(2,2)=strainl(2);
strtens(1,2)=0.5*strainl(3);
strtens(2,1)=0.5*strainl(3);
% if plane stress
if strcmp(mat.stressState,'PLANE_STRESS')
  strtens(3,3)=-mat.nu/(1-mat.nu)*(strainl(1)+strainl(2));
end

% compute principal strains
prstr          = eig(strtens);
prstr(prstr<0) = 0.;
eqvStrain      = norm(prstr);

% compute equivalent strain
% for i=1:3
%   if prstr(i) < 0.
%     prstr(i) = 0.;
%   end
% end
%eqeps=sqrt( prstr(1)^2. + prstr(2)^2. + prstr(3)^2. );



