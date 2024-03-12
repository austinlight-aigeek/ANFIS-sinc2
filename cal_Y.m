function Y = cal_Y(W, x1, x2, P, Q, R)

    if length(x1)~=length(x2) || length(x1)~=size(W,3)
        error('size of W and x1, x2 does not match!');
    end
    
    nRule = size(W,1)*size(W,2);

    ntrData = length(x1);
    A = zeros(ntrData, 3*nRule);
    
    for p = 1:ntrData
        Wp = W(:,:,p);
        A(p,:) = [mat2vector(Wp)'*x1(p) mat2vector(Wp)'*x2(p) mat2vector(Wp)'];
    end
    X = [mat2vector(P); mat2vector(Q); mat2vector(R)];
    Y = A*X;
end