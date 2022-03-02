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
        Li.*=0.0
        idx = A.colptr[i] : A.colptr[i + 1] - 1
        idx2 = L.colptr[i] : L.colptr[i + 1] - 1

        a.=0
        for j in idx
            #k+=1
            a[A.rowval[j]] = A.nzval[j]
        end

        @inbounds for j in idx2

            idx3 = L.colptr[L.rowval[j]] : L.colptr[L.rowval[j] + 1] - 1
            m = length(idx3)
            #Lj = view(L.nzval,idx3)
            #Lj.=0
            #Li1.=0
            Rdest = CartesianIndices((1:length(idx3)))
            Rsrc = CartesianIndices((idx3,),)
            #copyto!(Lj,Rdest,L.nzval,Rsrc)
            copy!(view(Lj,1:length(idx3)),view(L.nzval,idx3))

            #Lj = view(L,1:L.rowval[j],L.rowval[j])
            #println(i,"_",j)
            #println("---")
            #sd = view(L.rowval,idx3)
            #Rsrc = CartesianIndices((sd,))
            for (k,v) in enumerate(idx3)
                Li1[k] = Li[L.rowval[v]]
            end
            #Li1[1:length(idx3)] .= Li[view(L.rowval,idx3)]
            #copy!(view(Li1,1:length(idx3)),view(Li,sd))
            #copyto!(Li1,Rdest,L.nzval,Rsrc)
            sd = L.rowval[j]
            sdf = a[L.rowval[j]]
            hcho1(i,L,Li,Lj,sdf,sd,idx3,idx2,m);
        end
    end
    return nothing
end

function hcho1(i::Int64,L::SparseMatrixCSC{Float64, Int64},Li::Vector{Float64},
        Lj::Vector{Float64},val::Float64,j::Int64,idx3::UnitRange{Int64},idx2,m::Int64)
        s=0
        Li.=0
        Lj.=0
        Li[L.rowval[idx3]] = L.nzval[idx3]
        Lj[L.rowval[idx2]] = L.nzval[idx2]
        u = intersect(L.rowval[idx3],L.rowval[idx2])
        for k in 1:j
            #@inbounds s+=L.nzval[k[1]]*Li[L.rowval[k[2]]]
            @inbounds s+=Li[k]*Lj[k]
        end
        s=0
        for k in idx3
            if L.rowval[k] in L.rowval[idx2]
            #@inbounds s+=L.nzval[k[1]]*Li[L.rowval[k[2]]]
                @inbounds s+=L.nzval[k]*L.nzval[idx2[L.rowval[idx2].==L.rowval[k]]]
            end
        end
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
        Li[j] = L[j,i]
        return nothing
end
