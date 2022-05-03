function [CTO] = DP_CTO(matl,sigma,lambda,iter)

    if iter == 1
        CTO = matl.Ce;
    else
        [F,G] = DP_fg(sigma,matl,lambda);
        
        CTO   = matl.Ce - (matl.Ce*G.dG*F.dF'*matl.Ce)/...
                (F.dF'*matl.Ce*G.dG + F.fq'*F.ql);
    end

    return

end
