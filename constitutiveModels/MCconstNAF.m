function [sig,epsE,Dalg,yield_pos] = MCconstNAF(epsTr,E,v,phi,psi,c)

tol=1e-12; 
Ce=(-ones(3)*v+(1+v)*eye(3))/E; 
De=[inv(Ce) zeros(3); zeros(3) E/(2*(1+v))*eye(3)]; 
De3=De(1:3,1:3);

[t,epsVal]=eig([epsTr(1) epsTr(4)/2 epsTr(6)/2; epsTr(4)/2 epsTr(2) epsTr(5)/2; epsTr(6)/2 epsTr(5)/2 epsTr(3)]);
epsTr=[epsVal(1); epsVal(5); epsVal(9)];
[epsTr,sO]=sort(epsTr,'descend');
t=t(:,sO);
Q=[t(1)*t(1)   t(2)*t(2)   t(3)*t(3)   t(1)*t(2)           t(2)*t(3)           t(3)*t(1)           ;
   t(4)*t(4)   t(5)*t(5)   t(6)*t(6)   t(4)*t(5)           t(5)*t(6)           t(6)*t(4)           ;
   t(7)*t(7)   t(8)*t(8)   t(9)*t(9)   t(7)*t(8)           t(8)*t(9)           t(9)*t(7)           ;
   2*t(1)*t(4) 2*t(2)*t(5) 2*t(3)*t(6) t(1)*t(5)+t(4)*t(2) t(2)*t(6)+t(5)*t(3) t(3)*t(4)+t(6)*t(1) ;
   2*t(4)*t(7) 2*t(5)*t(8) 2*t(6)*t(9) t(4)*t(8)+t(7)*t(5) t(5)*t(9)+t(8)*t(6) t(6)*t(7)+t(9)*t(4) ;
   2*t(7)*t(1) 2*t(8)*t(2) 2*t(9)*t(3) t(7)*t(2)+t(1)*t(8) t(8)*t(3)+t(2)*t(9) t(9)*t(1)+t(3)*t(7)];
sig=Ce\epsTr; epsE=epsTr; Dalg=De;
k=(1+sin(phi))/(1-sin(phi)); sigC=2*c*sqrt(k);
f=k*sig(1)-sig(3)-sigC;
yield_pos=zeros(1,4);
if f>tol
  m=(1+sin(psi))/(1-sin(psi));  
  sigA=sigC/(k-1)*ones(3,1);
  r1=[1 1 k].'; r2=[1 k k].'; rg1=[1 1 m].'; rg2=[1 m m].';
  df=[k 0 -1].'; dg=[m 0 -1].'; 
  rp=De3*dg/(dg.'*De3*df);
  t1=(rg1.'*Ce*(sig-sigA))/(rg1.'*Ce*r1); 
  t2=(rg2.'*Ce*(sig-sigA))/(rg2.'*Ce*r2);
  f12=(cross(rp,r1)).'*(sig-sigA); f13=(cross(rp,r2)).'*(sig-sigA);
  if t1>tol && t2>tol        %Apex
    sig=sigA; Dep=zeros(6);
    yield_pos(1)=1;
  elseif f12<tol && f13<tol  %Line 1
    sig=sigA+t1*r1;
    Dep=[r1*rg1'/(r1'*Ce*rg1) zeros(3); zeros(3) E/(2*(1+v))*eye(3)];
    yield_pos(2)=1;
  elseif f12>tol && f13>tol  %Line 2
    sig=sigA+t2*r2;
    Dep=[r2*rg2'/(r2'*Ce*rg2) zeros(3); zeros(3) E/(2*(1+v))*eye(3)];
    yield_pos(3)=1;
  else                       %Plane
    sig=sig-f*rp;
    Dep=De3-(De3*(dg*df.')*De3)/(df.'*De3*dg);
    Dep=[Dep zeros(3); zeros(3) E/(2*(1+v))*eye(3)];
    yield_pos(4)=1;
  end
  epsE=Ce*sig; Dalg=Dep;
  T=zeros(3); sigTr=Ce\epsTr;
  if abs((sigTr(1)-sigTr(2)))>1e-3; T(1,1)=(sig(1)-sig(2))/(sigTr(1)-sigTr(2)); end
  if abs((sigTr(2)-sigTr(3)))>1e-3; T(2,2)=(sig(2)-sig(3))/(sigTr(2)-sigTr(3)); end
  if abs((sigTr(1)-sigTr(3)))>1e-3; T(3,3)=(sig(1)-sig(3))/(sigTr(1)-sigTr(3)); end
  Dalg(4:6,4:6)=T*Dalg(4:6,4:6);
  Dalg=Q.'*Dalg*Q;
end
sig=Q.'*[sig; zeros(3,1)]; epsE=Q\[epsE; zeros(3,1)];