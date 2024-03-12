function [SW, SWF] = cal_SWF(W, F)

    P = size(W,3);
    SW = zeros(1,P);
    SWF = SW;
    
    for p = 1:P
        SW(p) = sum(sum(W(:,:,p)));
        SWF(p) = sum(sum(W(:,:,p).*F(:,:,p)));
    end
    
end
