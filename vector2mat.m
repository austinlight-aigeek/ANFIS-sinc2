function M = vector2mat(m, ncol)
    nrow = length(m)/ncol;
    M = zeros(nrow, ncol);
    for i = 1:length(m)
        if mod(i,ncol)~=0
            M(ceil(i/ncol), mod(i,ncol)) = m(i);
        else
            M(ceil(i/ncol), ncol) = m(i);
        end
    end
end