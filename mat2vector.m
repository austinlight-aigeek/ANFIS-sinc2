function v = mat2vector(V)
    I = size(V,1);
    J = size(V,2);
    n = I*J;
    v = zeros(n,1);
    for i = 1:I
        v((i-1)*J+1:i*J) = V(i,:);
    end
end