function [W, W_bar] = cal_W(alpha1, alpha2, x1, x2, muAx)
    P = length(x1);
    
    W = zeros(4,4,P);
    W_bar = W;
    for p = 1:P
        for k = 1:4
            for m = 1:4
                W(k,m,p) = muAx(x1(p), alpha1(1,k),alpha1(2,k),alpha1(3,k))*...
                    muAx(x2(p), alpha2(1,m),alpha2(2,m),alpha2(3,m));
            end
        end
        W_bar(:,:,p) = W(:,:,p)/sum(sum(W(:,:,p)));
    end
end