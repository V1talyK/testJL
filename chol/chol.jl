function hand_cholS(A)
    L = sparse(LowerTriangular(A))
    L.nzval.=0.0
    n = size(A,1)
    L[1,1] = sqrt(A[1,1])
    Lr = zeros(Int64,Int64(n*n/2))
    Lc = zeros(Int64,0)
    Lv = zeros(Float64,0)

    Lr2 = zeros(Int64,0)
    Lc2 = zeros(Int64,0)
    Lv2 = zeros(Float64,0)

    LL = Vector{Array{Float32,1}}(undef,n)
    LL2 = Vector{Array{Float32,1}}(undef,n)
    for i=1:n
        LL[i] = zeros(Float32, n)
        LL2[i] = zeros(Float32, n)
    end

    Ld = zeros(n)

    Lr[1] = 1
    push!(Lc,1)
    push!(Lv,sqrt(A.nzval[1]))

    Ld[1] = sqrt(A.nzval[1])
    kr = 1
    for i = A.colptr[1]+1:A.colptr[2]-1
        #L[A.rowval[i],1] = A[A.rowval[i],1]/L[1,1]
        L[A.rowval[i],1] = A.nzval[i]/Ld[1]
        kr+=1
        Lr[kr] = A.rowval[i]
        push!(Lv,A.nzval[i]/Lv[1])

        LL[A.rowval[i]][1] = A.nzval[i]/Lv[1]
    end
    push!(Lc,kr+1)

    for j = 2:n
        s = 0.0
        #s0 = 0.0
        # for k = 1:j-1
        #     s+=L[j,k]^2
        #     #s+=v^2
        # end
        #s = 0
        # for k = Lc[1]:Lc[j]-1
        #     if Lr[k]==j
        #         s+=Lv[k]^2
        #     end
        # end
        s = sum(LL[j][1:j].*LL[j][1:j])
        #s0=sum(Lv[j].^2)
        #println(j," ",s," ",s0)
        for i = A.colptr[j]:A.colptr[j+1]-1
            if A.rowval[i]==j
                Ld[j] = sqrt(A.nzval[i] - s)
                kr+=1
                Lr[kr] = j
                #push!(Lr,j)
                push!(Lv,sqrt(A.nzval[i] - s))
                LL[A.rowval[i]][j] = sqrt(A.nzval[i] - s)
            end
        end
        # for i = A.colptr[j]:A.colptr[j+1]-1
        #     if A.rowval[i]>j
        #         s0 = 0.0
        #         for k = 1:j-1
        #             s0+=L[A.rowval[i],k]*L[j,k]
        #         end
        #         L[A.rowval[i],j] = (A.nzval[i]-s0)/L[j,j]
        #     end
        # end
        for ic = A.colptr[j]:A.colptr[j+1]-1
            if A.rowval[ic]>j
                LL[A.rowval[ic]][j] = A.nzval[ic]
            end
        end
        for i = j+1:n#A.colptr[j]:A.colptr[j+1]-1
            s = LinearAlgebra.BLAS.dot(j-1,LL[i],1,LL[j],1)
            s = 0.0
            @inbounds for k=1:j-1
                s+=LL[i][k]*LL[j][k]
            end
            s=LL2[i][j]+ LL[i][j-1]*LL[j][j-1]
            LL2[i][j] = s
            kr+=1
            Lr[kr] = i
            #push!(Lr,i)

            s1 = LL[i][j]
            LL[i][j] = 1/Ld[j]*(s1-s)
            push!(Lv,LL[i][j])
            #push!(Lv[i],1/Ld[j]*(A[i,j]-s))
        end
        push!(Lc,kr+1)
    end

    L[1:n+1:end].=Ld
    Lr = Lr[1:kr]
    return L, (Lr,Lc,Lv)
end

A.colptr

    A.rowval

@profiler hand_cholS(mA)
