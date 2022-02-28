A = sprand(10,10,0.2)
A = .-(A+A')
crc = CartesianIndex.(1:10,1:10)
A[crc].=-sum(A,dims=2).+1
A
b = zeros(10); b[1:5:10] .= -1
x = A\b

@time cholesky(A)

A.rowval
n = A.n
s = zeros(n)
L = zeros(n,n)
for i = 1:n
    for j = 1:i
        s = 0
        for k = 1:j
            s += L[i,k]*L[j,k]
        end
        if i == j
            L[i,j] = sqrt(A[i,i] - s);
        else
            L[i,j] = (1.0 / L[j,j] * (A[i,j] - s));
        end
    end
end

y = L\b
x1 = L'\y
