function normal = getPrincipalDirection( stress )

% d = stress(1) - stress(2);
%
% if ( abs ( d ) < 1e-6 )
%   theta = 0.78539816339744830962;
% else
%   theta = 0.5 * atan2 ( stress(3), 0.5 * d );
% end

% S = zeros(3,3);
% S(1,1) = stress(1);
% S(2,2) = stress(2);
% S(1,2) = stress(3);
% S(2,1) = stress(3);
%
% [V,D]  = eig(S);
% [~,I]  = sort(diag(D));
% normal = V(1:2,I(3));

sigmaxx = stress(1);
sigmayy = stress(2);
sigmaxy = stress(3);

d0      = (sigmaxx-sigmayy)^2+4.*sigmaxy^2;
d       = sqrt(d0);
temp    = 0.5*(sigmaxx+sigmayy);

sigma1 = temp + 0.5*d;
sigma2 = temp - 0.5*d;
tiny=1e-16;
if ((abs(sigma1)>tiny)||(abs(sigma2)>tiny))
  if (abs(sigmaxx-sigma1)>tiny)
    t = - (sigmaxy/(sigmaxx-sigma1)) ;
    normal(2) = 1.0 / sqrt(t*t + 1.0) ;
    normal(1) = t*normal(2) ;
  else
    if (abs(sigmayy-sigma1)>tiny)
      t = - (sigmaxy/(sigmayy-sigma1)) ;
    elseif (abs(sigmaxy)>tiny)
      t = - ((sigmaxx-sigma1)/sigmaxy) ;
    else 
      error('error in computing principal direction')      
    end
    normal(1) = 1.0 / sqrt(t*t + 1.0) ;
    normal(2) = t*normal(1) ;
  end
end

theta = acos(normal(1));

if ( theta > pi/2 ) || ( pi < -pi/2 )
  normal = - normal;
end

%normal(1)=abs(normal(1))
