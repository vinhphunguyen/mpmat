% An example of the two-scale cohesive model.
% Bulk: elastic
% Crack: initially rigid exponential cohesive law.
% VP Nguyen
% University of Adelaide, Australia, August, 2014.

E  = 10;
nu = 0.;
ft = 1;
Gf = 1;

H  = 11;

uf = 2*Gf/ft;

cracked = 0;

ao = elasticityMatrix(E,nu,'PLANE_STRAIN');

nstep = 100;

epsilon  = 0; sigma=0; u = 0; epsilon0=0;
epsilonf = 1;
dEps     = 1e-4;

K  = [0 0; 0 0];

n  = [1 0; 0 0; 0 1];
Ao = n'*ao*n;


stress  = [];
strain  = [];
strain0 = [];

jump   = [];
trac   = [];

while (abs(epsilon) < epsilonf) 
  if (cracked==0)
    sigma   = ao * [epsilon;0;0];
    yield = sigma(1) - ft;
    if (yield>=0) 
      cracked=1; 
      sigma(1)
      K(1,1) = -ft^2/Gf;
      %dEps = -dEps; % only for H=11 when snapback is to be modeled.
    end
    epsilon0 = epsilon;    
  else
    C  = K + (1/H)*Ao;
    DeltaU = inv(C)*n'*(ao*[dEps;0;0])
    if isnan(DeltaU(1)) 
      K(1,1)
      error('dd')
    end
    if (DeltaU(1)<0)
      dEps = -dEps;
      DeltaU = inv(C)*n'*(ao*[dEps;0;0])
      %error('negative jump')
    end
    epsilon0 = epsilon0 + [dEps;0;0]-(1/H)*n*DeltaU;
    epsilon0 = epsilon0(1);
    DeltaS = ao*([dEps;0;0]-(1/H)*n*DeltaU);
    u = u + DeltaU(1);
    t = ft*exp(-(ft/Gf)*u);
    K(1,1)=-ft^2/Gf*exp(-(ft/Gf)*u);
    sigma  = sigma + DeltaS;
    
    trac  = [trac;t];
    jump  = [jump;u];
  end
  
  stress = [stress;sigma(1)];
  strain = [strain;epsilon];
  strain0 = [strain0;epsilon0];
  % increasing load 
  epsilon = epsilon + dEps; 
end

save('tryH10Exp.mat','stress','strain','jump','trac');

h1=load('tryH1Exp.mat');
h2=load('tryH2Exp.mat');
h4=load('tryH4Exp.mat');
h8=load('tryH10Exp.mat');

figure
set(gca,'FontSize',14)
hold on
plot(h1.strain,h1.stress,'b-','LineWidth',1.6);
plot(h2.strain,h2.stress,'r-','LineWidth',1.6);
plot(h4.strain,h4.stress,'black-','LineWidth',1.6);
plot(h8.strain,h8.stress,'black--','LineWidth',1.6);
%plot(strain0,stress,'black-','LineWidth',1.6);
%plot(jump,trac,'r-','LineWidth',1.6);
xlabel('strain')
ylabel('stress')
legend('H=1','H=2','H=4','H=11')

%%

figure
set(gca,'FontSize',14)
hold on
plot(h1.jump,h1.trac,'b.-','LineWidth',1.6);
plot(h2.jump,h2.trac,'r--','LineWidth',1.6);
plot(h4.jump,h4.trac,'cy*','LineWidth',0.3);
plot(h8.jump,h8.trac,'black*','LineWidth',0.3);
%plot(strain0,stress,'black-','LineWidth',1.6);
%plot(jump,trac,'r-','LineWidth',1.6);
xlabel('jump')
ylabel('traction')
legend('H=1','H=2','H=4','H=11')


