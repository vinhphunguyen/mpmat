E  = 10;
nu = 0.;
ft = 1;
Gf = 1;

H  = 4;

uf = 2*Gf/ft;

cracked = 0;

ao = elasticityMatrix(E,nu,'PLANE_STRAIN');

nstep = 100;

epsilon  = 0; sigma=0; u = 0; epsilon0=0;
epsilonf = 1;
dEps     = 1e-4;

Kn = -ft/uf;
K  = [Kn 0; 0 0];

n  = [1 0; 0 0; 0 1];
Ao = n'*ao*n;
C  = K + (1/H)*Ao;

stress  = [];
strain  = [];
strain0 = [];

jump   = [];
trac   = [];

while (epsilon < epsilonf) 
  if (cracked==0)
    sigma   = ao * [epsilon;0;0];
    yield = sigma(1) - ft;
    if (yield>=0) 
      cracked=1; 
      sigma(1)
    end
    epsilon0 = epsilon;    
  else
    DeltaU = inv(C)*n'*(ao*[dEps;0;0]);
    epsilon0 = epsilon0 + [dEps;0;0]-(1/H)*n*DeltaU;
    epsilon0 = epsilon0(1);
    DeltaS = ao*([dEps;0;0]-(1/H)*n*DeltaU);
    u = u + DeltaU(1);
    t = ft*(1-u/uf);
    if (u>uf) t=0; C=(1/H)*Ao; end
    sigma  = sigma + DeltaS;
    
    trac  = [trac;t];
    jump  = [jump;u];
  end
  
  stress = [stress;sigma(1)];
  strain = [strain;epsilon];
  strain0 = [strain0;epsilon0];
  epsilon = epsilon + dEps;    
end

%save('tryH4.mat','stress','strain','jump','trac');

h1=load('tryH1.mat');
h2=load('tryH2.mat');
h4=load('tryH4.mat');

figure
set(gca,'FontSize',14)
hold on
plot(h1.strain,h1.stress,'b-','LineWidth',1.6);
plot(h2.strain,h2.stress,'r-','LineWidth',1.6);
plot(h4.strain,h4.stress,'black-','LineWidth',1.6);
%plot(strain0,stress,'black-','LineWidth',1.6);
%plot(jump,trac,'r-','LineWidth',1.6);
xlabel('strain')
ylabel('stress')
legend('H=1','H=2','H=4')

%%

figure
set(gca,'FontSize',14)
hold on
plot(h1.jump,h1.trac,'b.-','LineWidth',1.6);
plot(h2.jump,h2.trac,'r--','LineWidth',1.6);
plot(h4.jump,h4.trac,'cy*','LineWidth',0.3);
%plot(strain0,stress,'black-','LineWidth',1.6);
%plot(jump,trac,'r-','LineWidth',1.6);
xlabel('jump')
ylabel('traction')
legend('H=1','H=2','H=4')


