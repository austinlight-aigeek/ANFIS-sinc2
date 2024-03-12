function [delta_alpha_b1, delta_alpha_b2] = update_param_b(alpha,x1,x2,eta,E,W,F,muAx, DmuAx_b)
    
    [SW, SWF] = cal_SWF(W,F);
    alpha1 = alpha(:,:,1);
    alpha2 = alpha(:,:,2);
    P = length(SW);
    
    delta_alpha_b1 = zeros(1,4);
    delta_alpha_b2 = zeros(1,4);
    
    for p = 1:P
        Ep = E(p);
        Fp = F(:,:,p);
        SWp = SW(p);
        SWFp = SWF(p);
        x1p = x1(p);
        x2p = x2(p);
        D_alpha_b1 = zeros(1,4);
        D_alpha_b2 = zeros(1,4);
        for k = 1:4
            a1k = alpha1(1,k); a2k = alpha2(1,k);
            b1k = alpha1(2,k); b2k = alpha2(2,k);
            c1k = alpha1(3,k); c2k = alpha2(3,k);
            
            Db1WFp = 0; Db1Wp = 0;
            Db2WFp = 0; Db2Wp = 0;
            for j = 1:4
                a2j = alpha2(1,j); a1j = alpha1(1,j);
                b2j = alpha2(2,j); b1j = alpha1(2,j);
                c2j = alpha2(3,j); c1j = alpha1(3,j);
                
                Db1WFp = Db1WFp + muAx(x2p,a2j,b2j,c2j)*DmuAx_b(x1p,a1k,b1k,c1k)*Fp(k,j);
                Db1Wp = Db1Wp + muAx(x2p,a2j,b2j,c2j)*DmuAx_b(x1p,a1k,b1k,c1k);
                
                Db2WFp = Db2WFp + muAx(x1p,a1j,b1j,c1j)*DmuAx_b(x2p,a2k,b2k,c2k)*Fp(j,k);
                Db2Wp = Db2Wp + muAx(x1p,a1j,b1j,c1j)*DmuAx_b(x2p,a2k,b2k,c2k);
            end
            D_alpha_b1(k) = Db1WFp*SWp - Db1Wp*SWFp;
            D_alpha_b2(k) = Db2WFp*SWp - Db2Wp*SWFp;
            
        end
        
        delta_alpha_b1 = delta_alpha_b1 + Ep/SWp^2*D_alpha_b1;
        delta_alpha_b2 = delta_alpha_b2 + Ep/SWp^2*D_alpha_b2;
    end
    delta_alpha_b1 = eta*delta_alpha_b1;
    delta_alpha_b2 = eta*delta_alpha_b2;
    
    for i = 1:4
        if (delta_alpha_b1(i) - floor(delta_alpha_b1(i)))>0.5
            delta_alpha_b1(i) = ceil(delta_alpha_b1(i));
        elseif (delta_alpha_b1(i) - floor(delta_alpha_b1(i)))<0.5
            delta_alpha_b1(i) = floor(delta_alpha_b1(i));
        end

        if (delta_alpha_b2(i) - floor(delta_alpha_b2(i)))>0.5
            delta_alpha_b2(i) = ceil(delta_alpha_b2(i));
        elseif (delta_alpha_b2(i) - floor(delta_alpha_b2(i)))<0.5
            delta_alpha_b2(i) = floor(delta_alpha_b2(i));
        end
    end

end