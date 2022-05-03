function [sigma, alpha, epsilonp] = updateVonMisesAxisSym (epsilon,epsilonp0,alpha0,mu,kappa,k1,k0)

% update stress for a linear hardening von Mises material for axisymmetric
% problems.
% 
% Strain/stress stored as: [epsilon_xx epsilon_yy epsilon_xy epsilon_tt]
% i.e. the hoop component is the last entry.
%
% VP Nguyen, August 2015. 
% The Univerisity of Adelaide, Australia.

eye3   = [1 1 0 1]';
eye3x3 = [1 1 0 1; 
          1 1 0 1; 
          0 0 0 0;
          1 1 0 1];
I_dev  = eye(4) - (1/3)*eye3x3;
I      = [1 0 0 0; 
          0 1 0 0; 
          0 0 0.5 0; %0.5 to make engineering strain to physical one
          0 0 0 1]; 
Iinv   = [1 0 0 0; 0 1 0 0; 0 0 2 0; 0 0 0 1];

% Compute trial stress
epsilon_dev  = I_dev*epsilon;
s_trial      = 2*mu*I*(epsilon_dev-epsilonp0);
norm_s_trial = sqrt(s_trial(1)^2 + s_trial(2)^2 + s_trial(4)^2 + 2*s_trial(3)^2);
sigma_trial  = kappa*sum(epsilon(1:2)+epsilon(4))*eye3 + s_trial;

% Check yield condition
f_trial = norm_s_trial - (k1*alpha0 + k0);

if f_trial <= 0 % elastic step       
    alpha    = alpha0;
    epsilonp = epsilonp0;
    sigma    = sigma_trial;
else % plastic step
    normal = s_trial/norm_s_trial;
    lambda = (norm_s_trial - k1*alpha0 - k0)/(2*mu + k1);
    alpha = alpha0 + lambda;
    % Update plastic strain and stress
    epsilonp = epsilonp0 + lambda*Iinv*normal; 
    sigma = sigma_trial + s_trial - 2*mu*lambda*normal;
end
