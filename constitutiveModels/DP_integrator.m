function [sigma,epsp,lambda_n,iter] = DP_integrator(matl,eps,epsp,lambda)
    
    iter = 0;
    TOL = 10e-10;
    
    epsp_tr = epsp;
    dlambda = 0;
    
    for iter = 1:200
        sigma_tr = matl.Ce*(eps-epsp_tr)';
        [F,G]    = DP_fg(sigma_tr,matl,lambda);
        
        % Check if elastic
        if (F.F <= TOL)
            sigma    = sigma_tr;
            epsp     = epsp_tr;
            lambda_n = lambda;
            
            % disp(['Inner: ' num2str(iter)])
            
            return;
        end
            
        % Inelastic Material Response
        deltalambda = F.F/(F.dF'*matl.Ce*G.dG + F.fq'*F.ql);
        
        dlambda     = dlambda + deltalambda;
        lambda      = lambda   + deltalambda;
        epsp_tr     = epsp_tr  + (deltalambda*G.dG)';
    end
    
    disp(['Warning'])
    
    sigma    = sigma_tr;
    epsp     = epsp_tr;
    lambda_n = lambda;
    
    return
end
