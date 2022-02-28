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
for i = 1:n
    for j = 1:i
        s = 0
        k = 1:j
        s = dot(L[i,k],L[j,k])

        if i == j
            L[i,j] = sqrt(A[i,i] - s);
        else
            L[i,j] = (1.0 / L[j,j] * (A[i,j] - s));
        end
    end
end
A.rowval
A.colptr

function hcho(L::Matrix{Float64},A,Li::Vector{Float64},Lj::Vector{Float64},n,a)
    fl = falses(n)
    for i = 1:n
        Li.*=0.0
        idx = A.colptr[i] : A.colptr[i + 1] - 1

        a.=0
        for j in idx
            #k+=1
            a[A.rowval[j]] = A.nzval[j]
        end

        @inbounds for j = 1:i
            #println(i," ",j)
            fl.=false
            hcho1(i,L,Li,Lj,a[j],j,fl)
        end
    end
end

function hcho1(i,L,Li,Lj,val,j,fl)
        #copyto!(view(Lj,1:j),view(L,1:j,j))
        #copyto!(Lj,Rdest,Lj,Rsrc)
        #LinearAlgebra.BLAS.blascopy!(n,L[j,:],1,Lj,1)
        # s=0
        # @inbounds for k=1:j
        #     s+=Li[k]*L[j,k]
        # end
        #Lj = view(L,1:j,j)
        #s = LinearAlgebra.BLAS.dot(j,Li,1,Lj,1)
        Lj = view(L,1:j,j)
        s = LinearAlgebra.BLAS.dot(j,Li,1,Lj,1)
        fl[1:j] = view(L,1:j,j).!=0
        fl[1:j] = Li[1:j].!=0
        s = 0
        if any(fl)
            sd = view(L,:,j)[fl]
            df = Li[fl]
            s = LinearAlgebra.BLAS.dot(df,sd)
        end
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
