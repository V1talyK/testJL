n = mA.n
L = zeros(n,n)
for i = 1:n
    L[i,j] = sqrt(A[i,i])
    L[i,1] = A[i,1]/L[i,i]
    for j = 1:n

    end
end

mA.rowval
s = zeros(n)
for i = 1:1000
    for j = 1:i
        s = 0
        k = 1:j
        s = dot(L[i,k],L[j,k])

        if i == j
            L[i,j] = sqrt(mA[i,i] - s);
        else
            L[i,j] = (1.0 / L[j,j] * (mA[i,j] - s));
        end
    end
end
