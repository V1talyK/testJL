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
A.rowval
A.colptr

function hcho(L,A,Li::Vector{Float64},Lj::Vector{Float64},n,a)
    Li1 = zeros(n)
    for i = 1:n
        Li=L[:,i]
        idx = A.colptr[i] : A.colptr[i + 1] - 1
        idx2 = L.colptr[i] : L.colptr[i + 1] - 1

        a.=0
        for j in idx
            #k+=1
            a[A.rowval[j]] = A.nzval[j]
        end

        @inbounds for j in idx2
            Li = L[:,i]
            idx3 = L.colptr[L.rowval[j]] : L.colptr[L.rowval[j] + 1] - 1
            m = length(idx3)
            Rdest = CartesianIndices((1:length(idx3)))
            Rsrc = CartesianIndices((idx3,),)
            sd = L.rowval[j]
            sdf = a[L.rowval[j]]
            hcho1(i,L,Li,Lj,sdf,sd,idx3,idx2,m);
        end
    end
    return nothing
end

function hcho1(i::Int64,L::SparseMatrixCSC{Float64, Int64},Li,
        Lj,val::Float64,j::Int64,idx3::UnitRange{Int64},idx2,m::Int64)
        s=0
        #Li.=0
        Lj = L[:,j]
        #Li[view(L.rowval,idx3)] = view(L.nzval,idx3)
        #Lj[view(L.rowval,idx2)] = view(L.nzval,idx2)
        # u = intersect(L.rowval[idx3],L.rowval[idx2])
        for k in 1:j
            #@inbounds s+=L.nzval[k[1]]*Li[L.rowval[k[2]]]
            #@inbounds s+=Li[k]*Lj[k]
            #@inbounds s+=Li[k]*Lj[k]
        end
        s = dot(Li[1:j],Lj[1:j])
        # s=0
        # for k = idx3
        #     for k1 = idx2
        #     @inbounds if isequal(L.rowval[k], L.rowval[k1])
        #     #@inbounds s+=L.nzval[k[1]]*Li[L.rowval[k[2]]]
        #         s+=L.nzval[k]*L.nzval[k1]
        #     end
        #     end
        # end
        # for k in 1:j
        #     @inbounds s+=L[k,i]*L[k,j]
        # end
        # println(L.nzval[idx3])
        # println(L.nzval[idx2])
        # println(Vector(L[1:j,i]).*Vector(L[1:j,j]))
        # #println()
        # println("__________")
        #s = LinearAlgebra.BLAS.dot(m,Lj,1,Li1,1)
        #a1 = Lj.*cLi
        #s = reduce(+,a1)
        #s = LinearAlgebra.BLAS.dot(j,Li,1,Lj,1)
        s = val - s

        if i == j
            L[j,i] = sqrt(s);
        else
            L[j,i] = s / L[j,j];
        end
        #Li[j] = L[j,i]
        return nothing
end

function chol_fr(A,L)
    for i = 1:A.n
        L[i,i] = sqrt(A[i, i])
        for j = i+1:n
            L[j,i] = L[j,i]/L[i,i]
        end
        for k=i+1:n
            for j = k:n
                L[j,k] = A[j,k] - L[j,i]*L[k,i]
            end
        end
    end
end

chol_fr(A,L1)

y = L1\b;
x2 = L1'\y;
sum(abs,x0.-x2)
