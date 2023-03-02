function hand_cholS(A)
    L = spzeros(size(A))
    n = size(A,1)
    L[1,1] = sqrt(A[1,1])
    for i = A.colptr[1]+1:A.colptr[2]-1
        L[A.rowval[i],1] = A[A.rowval[i],1]/L[1,1]
    end

    for j = 2:n
        s = 0.0
        for k=1:j-1
            s+=L[j,k]^2
        end
        L[j,j] = sqrt(A[j,j] - s)
        for i = j+1:n#A.colptr[j]:A.colptr[j+1]-1
            s = 0.0
            for k = 1:j-1
                s += L[i,k]*L[j,k]
            end
            L[i,j] = 1/L[j,j]*(A[i,j]-s)
        end
    end
    return L
end

A.colptr

    A.rowval

@profiler hand_cholS(mA)
