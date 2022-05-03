function [matl] = DP_init(E,nu,c,phi,psi,Hci,Hpi,Hsi)
I = [1 1 1 0 0 0]';     

matl.E     = E;     
matl.nu    = nu;         
matl.c0    = c;      
matl.phi0  = phi;
matl.psi0  = psi;

matl.Hci   = Hci; % -0.0e-2;
matl.Hpi   = Hpi; % -5.0e-2;
matl.Hsi   = Hsi; % -0.0e-1;

matl.phi   = matl.phi0*pi()/180;
matl.psi   = matl.psi0*pi()/180;

% 2nd order tensor of elastic constants. Isotropic linear elastic.
temp    = matl.E/((1+matl.nu)*(1-2*matl.nu));
nu      = matl.nu;
matl.Ce = temp * ...
          [1-nu nu nu 0 0 0;
           nu 1-nu nu 0 0 0;
           nu nu 1-nu 0 0 0;
           0 0 0 0.5-nu 0 0;
           0 0 0 0 0.5-nu 0;
           0 0 0 0 0 0.5-nu];
       
matl.Hc = matl.Hci * norm(matl.Ce,2);        
matl.Hp = matl.Hpi * norm(matl.Ce,2);
matl.Hs = matl.Hsi * norm(matl.Ce,2);

matl.h = 1;
matl.H = 2;
matl.f = matl.h/matl.H;

return
end
