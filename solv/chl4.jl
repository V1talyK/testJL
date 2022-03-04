function hcho3(L,A,Li::Vector{Float64},Lj::Vector{Float64},n,a)
    AA = A.*A'
    for i = 1:n
        Li=L[:,i]
        idx = A.colptr[i] : A.colptr[i + 1] - 1
        idx2 = L.colptr[i] : L.colptr[i + 1] - 1
        a.=0
        for j in idx
            a[A.rowval[j]] = A.nzval[j]
        end

        @inbounds for j in idx2
            Li = L[:,i]
            idx3 = L.colptr[L.rowval[j]] : L.colptr[L.rowval[j] + 1] - 1
            m = length(idx3)
            sd = L.rowval[j]
            sdf = a[L.rowval[j]]
            hcho2(i,L,Li,Lj,sdf,sd,idx3,idx2,m, AA);
        end
    end
    return nothing
end

function hcho2(i::Int64,L::SparseMatrixCSC{Float64, Int64},Li,
        Lj,val::Float64,j::Int64,idx3::UnitRange{Int64},idx2,m::Int64, AA)
        Lj = L[:,j]
        s = sum(AA[1:j,i])
        s = val - s

        if i == j
            L[j,i] = sqrt(s);
        else
            L[j,i] = s / L[j,j];
        end
        #Li[j] = L[j,i]
        return nothing
end
