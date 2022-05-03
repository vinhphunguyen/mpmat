function [F,G,p,q,alp1,k,c] = DP_fg(sigma,matl,lambda)
% sigma -> s
I  = [1 1 1 0 0 0]';
s1 = sigma(1);          % sig_11
s2 = sigma(2);          % sig_22
s3 = sigma(3);          % sig_33
s4 = sigma(4);          % sig_12
s5 = sigma(5);          % sig_23  
s6 = sigma(6);          % sig_31

I1  = (s1+s2+s3)/3;
J2  = ((s1-s2)^2+(s2-s3)^2+(s3-s1)^2+6*(s4^2+s5^2+s6^2))/6;
A   = sqrt(9*J2);

stt = +6*s1-3*s2-3*s3;
tst = -3*s1+6*s2-3*s3;
tts = -3*s1-3*s2+6*s3;

% Material parameters as a function of the accumulative plastic strain
Kc  = matl.Hc;
Kp  = matl.Hp;
Ks  = matl.Hs;
c   = matl.c0  + Kc*lambda;
phi = matl.phi + Kp*lambda;
psi = matl.psi + Ks*lambda;

sph = sin(phi); cph = cos(phi);
ssh = sin(psi); csh = cos(psi);
sq3 = sqrt(3);

% Three-dimensional
alp1   = 6*sph/(sq3*(3-sph));
k      = 6*cph/(sq3*(3-sph));
alp2   = 6*sin(psi)/(sq3*(3-sin(psi)));

% Yield function F
F.F   = alp1*I1 + A/3 - k*c;

F.dF  = [stt   ;...
         tst   ;...
         tts   ;...
         18*s4 ;...
         18*s5 ;...
         18*s6 ]./(6*A) + alp1*I/3;

F.fq  = [I1; -1];
F.ql  = [2*sq3*Kp*cph/(3-sph) + 2*sq3*Kp*cph*sph/((3-sph)^2); 
         2*sq3*matl.c0*cph/(3-sph) - 2*sq3*Kp*c*sph/(3-sph) + ...
         2*sq3*Kp*c*cph*cph/((3-sph)^2)];

% Plastic potential function G
G.G   = alp2*I1 + A/3;

G.dG  = [stt   ;...
         tst   ;...
         tts   ;...
         18*s4 ;...
         18*s5 ;...
         18*s6 ]./(6*A) + alp2*I/3;

G.d2G = [stt*stt   stt*tst   stt*tts   stt*18*s4 stt*18*s5 stt*18*s6 ;
         tst*stt   tst*tst   tst*tts   tst*18*s4 tst*18*s5 tst*18*s6 ;
         tts*stt   tts*tst   tts*tts   tts*18*s4 tts*18*s5 tts*18*s6 ;
         18*s4*stt 18*s4*tst 18*s4*tts 324*s4*s4 324*s4*s5 324*s4*s6 ;
         18*s5*stt 18*s5*tst 18*s5*tts 324*s5*s4 324*s5*s5 324*s5*s6 ;
         18*s6*stt 18*s6*tst 18*s6*tts 324*s6*s4 324*s6*s5 324*s6*s6 ;
        ]./(-12*A^3) + ...
        [+2 -1 -1  0  0  0 ;
         -1 +2 -1  0  0  0 ; 
         -1 -1 +2  0  0  0 ;
          0  0  0  6  0  0;
          0  0  0  0  6  0;
          0  0  0  0  0  6; ]./(2*A);

G.dGl = (2*sq3*Ks*csh/(3-ssh) + 2*sq3*Ks*csh*ssh/((3-ssh)^2))*I/3;      
      
p = I1;
q = A/sq3;
      
return
end
