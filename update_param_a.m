function [delta_alpha_a1, delta_alpha_a2] = update_param_a(alpha,x1,x2,eta,E,W,F,muAx, DmuAx_a)
    
    [SW, SWF] = cal_SWF(W,F);
    alpha1 = alpha(:,:,1);
    alpha2 = alpha(:,:,2);
    P = length(SW);
    
    delta_alpha_a1 = zeros(1,4);
    delta_alpha_a2 = zeros(1,4);
    
    for p = 1:P
        Ep = E(p);
        Fp = F(:,:,p);
        SWp = SW(p);
        SWFp = SWF(p);
        x1p = x1(p);
        x2p = x2(p);
        D_alpha_a1 = zeros(1,4);
        D_alpha_a2 = zeros(1,4);
        for k = 1:4
            a1k = alpha1(1,k); a2k = alpha2(1,k);
            b1k = alpha1(2,k); b2k = alpha2(2,k);
            c1k = alpha1(3,k); c2k = alpha2(3,k);
            
            Da1WFp = 0; Da1Wp = 0;
            Da2WFp = 0; Da2Wp = 0;
            for j = 1:4
                a2j = alpha2(1,j); a1j = alpha1(1,j);
                b2j = alpha2(2,j); b1j = alpha1(2,j);
                c2j = alpha2(3,j); c1j = alpha1(3,j);
                
                Da1WFp = Da1WFp + muAx(x2p,a2j,b2j,c2j)*DmuAx_a(x1p,a1k,b1k,c1k)*Fp(k,j);
                Da1Wp = Da1Wp + muAx(x2p,a2j,b2j,c2j)*DmuAx_a(x1p,a1k,b1k,c1k);
                
                Da2WFp = Da2WFp + muAx(x1p,a1j,b1j,c1j)*DmuAx_a(x2p,a2k,b2k,c2k)*Fp(j,k);
                Da2Wp = Da2Wp + muAx(x1p,a1j,b1j,c1j)*DmuAx_a(x2p,a2k,b2k,c2k);
            end
            D_alpha_a1(k) = Da1WFp*SWp - Da1Wp*SWFp;
            D_alpha_a2(k) = Da2WFp*SWp - Da2Wp*SWFp;
            
        end
        
        delta_alpha_a1 = delta_alpha_a1 + Ep/SWp^2*D_alpha_a1;
        delta_alpha_a2 = delta_alpha_a2 + Ep/SWp^2*D_alpha_a2;
    end
    delta_alpha_a1 = eta*delta_alpha_a1;
    delta_alpha_a2 = eta*delta_alpha_a2;

end