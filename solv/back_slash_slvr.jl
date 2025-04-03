function back_slash_slvr!(x,
                         A::NamedTuple{(:L, :U, :p, :x_temp), Tuple{SparseMatrixCSC{Float64, Int64},
                                                                SparseMatrixCSC{Float64, Int64},
                                                                Vector{Int64},
                                                                Matrix{Float64}}},
                        b)

    ldiv_cl!(x,A,b)
    return nothing
end

function back_slash_slvr!(x,A::SuiteSparse.CHOLMOD.Factor{Float64},b)
    x .= A\b;
    return nothing
end

function make_CL_in_julia(ACL, nth = 1)
    LL = sparse(ACL.L)
    UU = copy(LL')
    x_temp = zeros(LL.n, nth)
    return (L = LL, U = UU, p = ACL.p, x_temp)
end

# function updateCL!(CL, ACL)
#     @timeit to1 "6.1" CL.L.=sparse(ACL.L)
#     @timeit to1 "6.2" CL.U.=copy(CL.L')
# end

function ldiv_cl!(x,
                  CL::NamedTuple{(:L, :U, :p, :x_temp),
                      Tuple{SparseMatrixCSC{Float64, Int64},
                      SparseMatrixCSC{Float64, Int64},
                      Vector{Int64},
                      Matrix{Float64}}},
                  b)
    @inbounds x .= view(b,CL.p)
    @inbounds xp = view(x,CL.p)
    #x = x[CL.p]
    x_temp = view(CL.x_temp,:,Threads.threadid())
    #x_temp = CL.x_temp[:]
    forward_substit!(x_temp,CL.L, x)
    backward_substit!(x,CL.U, x_temp)
    copy!(x_temp, x)
    @inbounds x[CL.p] = x_temp
end

function forward_substit!(x, S, b)
    x .= 0;
     @fastmath @inbounds for col = 1:S.n
        tmp = S.colptr[col]
        xc = (b[col] + x[col])/S.nzval[tmp]
        idx = tmp+1 : S.colptr[col + 1] - 1
        #v1 = view(S.rowval,idx)
        for v in idx
             x[S.rowval[v]] -=  S.nzval[v] * xc
        end
        x[col] = xc
    end
end

function backward_substit!(x::Vector{Float64},
                           UU::SparseMatrixCSC{Float64, Int64},
                           b)
    x .= 0;
    tmp0 = UU.colptr[UU.n + 1]
    @inbounds for col = UU.n:-1:1
        tmp = UU.colptr[col]
        xc = (b[col] + x[col])/UU.nzval[tmp0-1]
         
        idx = tmp0-2 :-1 : tmp
        for v in idx
             x[UU.rowval[v]] -=  UU.nzval[v] * xc
        end
        tmp0 = tmp
        x[col] = xc
    end
end