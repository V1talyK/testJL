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
        ifl = findall(fl[1:j])
        if any(fl)
            sd = copy(view(L,ifl,j))
            df = copy(view(Li,ifl))
            s = LinearAlgebra.BLAS.dot(j,df,1,sd,1)
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
