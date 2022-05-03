function [bulk,shear] = computeBulkShear(E,nu)

shear  = E/2/(1+nu);                    % shear modulus
bulk   = E/3/(1-2*nu);
