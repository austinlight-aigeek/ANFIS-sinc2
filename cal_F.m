function F = cal_F(P,Q,R,x1,x2)
    nTrain = length(x1);
    F = zeros(4,4,nTrain);
    
    for p = 1:nTrain
        for i = 1:4
            for j = 1:4
                F(i,j,p) = P(i,j)*x1(p)+Q(i,j)*x2(p)+R(i,j);
            end
        end
    end
    
end