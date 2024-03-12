function [P, Q, R] = update_param2(W, x1, x2, yd)

	nRule = size(W,1)*size(W,2);

    P = length(x1);
    A = zeros(P, 3*nRule);
    
    for p = 1:P
        Wp = W(:,:,p);
        A(p,:) = [mat2vector(Wp)'*x1(p) mat2vector(Wp)'*x2(p) mat2vector(Wp)'];
    end
    lamda = 5;
    S = lamda*eye(3*nRule);
    X = zeros(3*nRule,1);
    
    for i = 0:P-1
        a = A(i+1,:);
        b = yd(i+1,:);
        
        S = S - (S*a'*a*S)/(1+a*S*a');
        X = X + S*a'*(b-a*X);
    end
    pqr = reshape(X, [], 3);
    P = vector2mat(pqr(:,1), size(W,2));
    Q = vector2mat(pqr(:,2), size(W,2));
    R = vector2mat(pqr(:,3), size(W,2));
end