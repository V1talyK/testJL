L2 = zeros(n,n)
@time for i = 1:n
    v = A.colptr[i]:A.colptr[i+1]-1
    for j = 1:i
        #if A.rowval[j]<=i
            s = 0
            k = 1:j
            #s = dot(view(L2,i,k),view(L2,j,k))
            if i == j
                L2[i,j] = sqrt(A[i,i] - s);
            else
                L2[i,j] = (1.0 / L2[j,j] * (A[i,j] - s));
            end
        #end
    end
end
y = L2\b
x2 = L2'\y

sum(abs,x2.-x)

B = LowerTriangular(A)
