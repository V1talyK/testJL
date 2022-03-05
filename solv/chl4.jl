function chol11(L,A,n)
    for i = 1:n
        idx1 = L.colptr[i] : L.colptr[i + 1] - 1
        for j in view(L.rowval,idx1)
            idx2 = L.colptr[j] : L.colptr[j + 1] - 1
            v = view(L.rowval,idx2)
            ia = indexin(view(L.rowval,idx1),v)
            ia = ia[.!isnothing.(ia)]
            s = 0
            k1=0
            for k in 1:length(ia)
                k1+=1
                s += L.nzval[idx1[ia[k]]]*L.nzval[idx2[k1]]
                #s += L[k,i]*L.nzval[idx2[k1]]
            end
            s = A[j,i] - s
            if i == j
                L[i,j] = sqrt(s);
            else
                L[j,i] = (s / L[j,j]);
            end
        end
    end
end


@time chol11(L,A,n)
@profiler  chol11(L,A,n)
@btime begin
    L.=0
    chol11($L,$A,$n)
end

y = L'\b;
    x1 = L\y;
    sum(abs,x0.-x1)
