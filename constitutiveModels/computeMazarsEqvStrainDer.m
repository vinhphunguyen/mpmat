function dEqvdStr = computeMazarsEqvStrainDer(strainl,mat)

poiss = mat.nu;

% form strain tensor
strtens=zeros(3,3);
strtens(1,1)=strainl(1);
strtens(2,2)=strainl(2);
strtens(1,2)=0.5*strainl(3);
strtens(2,1)=0.5*strainl(3);
% if plane stress
if strcmp(mat.stressState,'PLANE_STRESS')
  strtens(3,3)=-poiss/(1-poiss)*(strainl(1)+strainl(2));
end

% compute principal strains and principal directions (rotation matrix)
[rot,prstr]=eig(strtens);

% compute positive principal strains
for i=1:3
  if prstr(i,i) < 0.
    prstr(i,i) = 0.;
  end
end

% derivative of equivalent strain in principal strain space
den=sqrt( prstr(1,1)^2. + prstr(2,2)^2. + prstr(3,3)^2. );
if den == 0.
  den =1;
end

for i=1:3
  deqprin(i)= prstr(i,i) / den;
end

% form tensor of derivatives
dertens=zeros(3,3);
for i=1:3
  dertens(i,i)=deqprin(i);
end

% rotate derivatives of strain tensor back to global reference system
tensrot=rot*dertens*rot';

% extract vector of derivatives
deqdep(1)=tensrot(1,1);
deqdep(2)=tensrot(2,2);
deqdep(3)=0.5*tensrot(1,2);
dEqvdStr=deqdep';

%   fprintf(1,'1 deqdep  %12.8g %12.8g %12.8g \n',deqdep(1),deqdep(2),deqdep(3))

% delta=sqrt(strainl(1)^2-2*strainl(1)*strainl(2)+strainl(2)^2+strainl(3)^2);
% 
% if delta == 0.
%   delta =1;
% end
% 
% de1dex =0.5*(delta+strainl(1)-strainl(2))/delta;
% de1dey =0.5*(delta-strainl(1)+strainl(2))/delta;
% de1dgxy=0.5*strainl(3)/delta;
% 
% de2dex = 0.5*(delta-strainl(1)+strainl(2))/delta;
% de2dey = 0.5*(delta+strainl(1)-strainl(2))/delta;
% de2dgxy=-0.5*strainl(3)/delta;
% 
% de3dex = poiss/(poiss-1);
% de3dey = poiss/(poiss-1);
% de3dgxy= 0;
% 
% xdeqdep(1)=deqprin(1)*de1dex  + deqprin(2)*de2dex  + deqprin(3)*de3dex;
% xdeqdep(2)=deqprin(1)*de1dey  + deqprin(2)*de2dey  + deqprin(3)*de3dey;
% xdeqdep(3)=deqprin(1)*de1dgxy + deqprin(2)*de2dgxy + deqprin(3)*de3dgxy;
% xdeqdep=xdeqdep';


