function [delta_alpha_c1, delta_alpha_c2] = update_param_c(alpha,x1,x2,eta,E,W,F,muAx, DmuAx_c)
    
    [SW, SWF] = cal_SWF(W,F);
    alpha1 = alpha(:,:,1);
    alpha2 = alpha(:,:,2);
    P = length(SW);
    
    delta_alpha_c1 = zeros(1,4);
    delta_alpha_c2 = zeros(1,4);
    
    for p = 1:P
        Ep = E(p);
        Fp = F(:,:,p);
        SWp = SW(p);
        SWFp = SWF(p);
        x1p = x1(p);
        x2p = x2(p);
        D_alpha_c1 = zeros(1,4);
        D_alpha_c2 = zeros(1,4);
        for k = 1:4
            a1k = alpha1(1,k); a2k = alpha2(1,k);
            b1k = alpha1(2,k); b2k = alpha2(2,k);
            c1k = alpha1(3,k); c2k = alpha2(3,k);
            
            Dc1WFp = 0; Dc1Wp = 0;
            Dc2WFp = 0; Dc2Wp = 0;
            for j = 1:4
                a2j = alpha2(1,j); a1j = alpha1(1,j);
                b2j = alpha2(2,j); b1j = alpha1(2,j);
                c2j = alpha2(3,j); c1j = alpha1(3,j);
                
                Dc1WFp = Dc1WFp + muAx(x2p,a2j,b2j,c2j)*DmuAx_c(x1p,a1k,b1k,c1k)*Fp(k,j);
                Dc1Wp = Dc1Wp + muAx(x2p,a2j,b2j,c2j)*DmuAx_c(x1p,a1k,b1k,c1k);
                
                Dc2WFp = Dc2WFp + muAx(x1p,a1j,b1j,c1j)*DmuAx_c(x2p,a2k,b2k,c2k)*Fp(j,k);
                Dc2Wp = Dc2Wp + muAx(x1p,a1j,b1j,c1j)*DmuAx_c(x2p,a2k,b2k,c2k);
            end
            D_alpha_c1(k) = Dc1WFp*SWp - Dc1Wp*SWFp;
            D_alpha_c2(k) = Dc2WFp*SWp - Dc2Wp*SWFp;
            
        end
        
        delta_alpha_c1 = delta_alpha_c1 + Ep/SWp^2*D_alpha_c1;
        delta_alpha_c2 = delta_alpha_c2 + Ep/SWp^2*D_alpha_c2;
    end
    delta_alpha_c1 = eta*delta_alpha_c1;
    delta_alpha_c2 = eta*delta_alpha_c2;

end