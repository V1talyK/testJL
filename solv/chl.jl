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
    fl = falses(n)
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

        @inbounds @fastmath for j in idx2
            #println(i," ",j)
            fl.=false
            #r1 = L.rowval[j]
            idx3 = L.colptr[L.rowval[j]] : L.colptr[L.rowval[j] + 1] - 1
            #Lj = view(L.nzval,idx3)
            #Lj.=0
            #Li1.=0
            Rdest = CartesianIndices((1:length(idx3)))
            Rsrc = CartesianIndices((idx3,),)
            copyto!(Lj,Rdest,L.nzval,Rsrc)
            #copy!(view(Lj,1:length(idx3)),view(L.nzval,idx3))

            #Lj = view(L,1:L.rowval[j],L.rowval[j])
            #println(i,"_",j)
            #println("---")
            sd = view(L.rowval,idx3)
            #Rsrc = CartesianIndices((sd,))
            for (k,v) in enumerate(idx3)
                Li1[k] = Li[L.rowval[v]]
            end
            #Li1[1:length(idx3)] .= Li[view(L.rowval,idx3)]
            #copy!(view(Li1,1:length(idx3)),view(Li,sd))
            #copyto!(Li1,Rdest,L.nzval,Rsrc)
            hcho1(i,L,Li,Li1,Lj,a[L.rowval[j]],L.rowval[j],fl,sd)
        end
    end
end

function hcho1(i,L,Li,Li1,Lj,val,j,fl,idx3)
        #copyto!(view(Lj,1:j),view(L,1:j,j))
        #copyto!(Lj,Rdest,Lj,Rsrc)
        #LinearAlgebra.BLAS.blascopy!(n,L[j,:],1,Lj,1)
        # s=0
        # @inbounds for k=1:j
        #     s+=Li[k]*L[j,k]
        # end
        #s = LinearAlgebra.BLAS.dot(j,Li,1,Lj,1)
        #Lj = view(L,1:j,j)
        #s = sum(Lj.*view(Li,1:j))
        #@inbounds cLi = view(Li,idx3)
        # println(Lj[1:length(idx3)])
        # println(Li1[1:length(idx3)])
        # println(cLi)
        # println("---")
        #s = LinearAlgebra.dot(Lj,Li1)
        #s = LinearAlgebra.dot(Lj[1:length(idx3)],cLi)
        #s = LinearAlgebra.dot(Lj[1:length(idx3)],Li1[1:length(idx3)])
        s = LinearAlgebra.BLAS.dot(length(idx3),Lj,1,Li1,1)
        #a1 = Lj.*cLi
        #s = reduce(+,a1)
        #s = LinearAlgebra.BLAS.dot(j,Li,1,Lj,1)
        s = val - s

        if i == j
            L[j,i] = sqrt(s);
            Li[j] = L[j,i];
            #Lj[i] = L[i,j];
        else
            Li[j] = s / L[j,j];
            L[j,i] = Li[j]
            #Lj[i] = L[i,j]
        end
end
