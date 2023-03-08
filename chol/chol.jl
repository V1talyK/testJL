function hand_cholS(A)
    L = sparse(LowerTriangular(A))
    L.nzval.=0.0
    n = size(A,1)
    L[1,1] = sqrt(A[1,1])
    Lr = zeros(Int64,0)
    Lc = zeros(Int64,0)
    Lv = zeros(Float64,0)
    Ld = zeros(n)

    push!(Lr,1)
    push!(Lc,1)
    push!(Lv,sqrt(A.nzval[1]))
    Ld[1] = sqrt(A.nzval[1])

    for i = A.colptr[1]+1:A.colptr[2]-1
        #L[A.rowval[i],1] = A[A.rowval[i],1]/L[1,1]
        L[A.rowval[i],1] = A.nzval[i]/Ld[1]
        push!(Lr,A.rowval[i])
        push!(Lv,A.nzval[i]/Lv[1])
    end
    push!(Lc,length(Lr)+1)

    for j = 2:n
        s = 0.0
        #s0 = 0.0
        # for k = 1:j-1
        #     s+=L[j,k]^2
        #     #s+=v^2
        # end
        #s = 0
        for k = Lc[1]:Lc[j]-1
            if Lr[k]==j
                s+=Lv[k]^2
            end
        end
        #s0=sum(Lv[j].^2)
        #println(j," ",s," ",s0)
        for i = A.colptr[j]:A.colptr[j+1]-1
            if A.rowval[i]==j
                Ld[j] = sqrt(A.nzval[i] - s)
                push!(Lr,j)
                push!(Lv,sqrt(A.nzval[i] - s))
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
        for i = j+1:n#A.colptr[j]:A.colptr[j+1]-1
            s = 0.0
            # for k = 1:j-1
            #     s += L[i,k]*L[j,k]
            # end
            s = 0.0
            sj = 0
            si = 0
            for k1=1:j-1
                for k = Lc[k1]:Lc[k1+1]-1
                    if Lr[k]==j
                        sj = Lv[k]
                    end
                    if Lr[k]==i
                        si = Lv[k]
                    end
                end
                s+=si*sj
                si = 0
                sj = 0
            end

            #L[i,j] = 1/Ld[j]*(A[i,j]-s)
            push!(Lr,i)
            push!(Lv,1/Ld[j]*(A[i,j]-s))
            #push!(Lv[i],1/Ld[j]*(A[i,j]-s))
        end
        push!(Lc,length(Lr)+1)
    end

    L[1:n+1:end].=Ld
    return L, (Lr,Lc,Lv)
end

A.colptr

    A.rowval

@profiler hand_cholS(mA)
