function backward_substitution!(U, y::AbstractVector)
    n = size(U,1)
    @inbounds for col = n : -1 : 1

        # Substitutions
        for idx = U.data.colptr[col + 1] - 1 : -1 : U.data.colptr[col] + 1
            y[col] -= U.data.nzval[idx] * y[U.data.rowval[idx]]
        end

        # Final answer for y[col]
        y[col] /= U.data.nzval[U.data.colptr[col]]
    end

    y
end

function forward_substitution!(y, L, b::AbstractVector)
    y = copy(b)
    @inbounds for col = 1 : L.data.n - 1
        for idx = L.data.colptr[col] : L.data.colptr[col + 1] - 1
            y[L.data.rowval[idx]] -= L.data.nzval[idx] * y[col]
        end
    end

    y
end


function forward_substitution1!(L, y::AbstractVector)
    @inbounds for col = 1 : L.data.n - 1
        # for idx = L.data.colptr[col] : L.data.colptr[col + 1] - 1
        #     y[L.data.rowval[idx]] -= L.data.nzval[idx] * y[col]
        # end
        idx = L.data.colptr[col] : L.data.colptr[col + 1] - 1
        iy = view(L.data.rowval,idx)
        z = view(L.data.nzval,idx)
        #y[iy] .= view(y,iy) .- z.*y[col]
        nmn!(view(y,iy),z,-y[col])
    end

    y
end

function nmn!(a,b,c)
    axpy!(c, b, a)
end

@btime $x .= $ACL\$b

LL = LowerTriangular(sparse(ACL.L))
bp = copy(b[ACL.p])
UU = copy(transpose(LL))


@btime y = backward_substitution!(LL, bp)
y2 = copy(y)
y1 = copy(y)

@btime x1 = forward_substitution!(UU, y1)
@btime x2 .= forward_substitution1!(UU, y2)
@profiler forward_substitution1!(UU, y2)

sum(abs,x[CL.p].-x1)

@btime ldiv!($y,$LL,$bp)
@btime ldiv!($x1,$UU,$y)


function frw_sb0!(x, S, b)
    x .= 0;
    @inbounds for col = 1:S.n
        idx = S.colptr[col] : S.colptr[col + 1] - 1
        #println(idx," ",S.rowval[idx])
        @inbounds for (k,v) in enumerate(idx)
            #println(j)
            j = S.rowval[v]
            if j == col
                x[j] = (b[j] + x[j])/S.nzval[idx[k]]
            else
                x[j] -=  S.nzval[idx[k]] * x[col]
            end
        end
    end
    return x
end


function frw_sbt!(x, S, b)
    x .= 0;
    @inbounds for col = 1:S.n
        idx = S.colptr[col]+1 : S.colptr[col + 1] - 1
        #println(idx," ",S.rowval[idx])
        x[col] = (b[col] + x[col])/S.nzval[S.colptr[col]]
        v1 = view(S.rowval,idx)
        #for i in eachindex(idx,v1)
             #j = S.rowval[v]
        @tturbo x[v1] .= view(x,v1) .- view(S.nzval,idx) .* x[col]
        #fg.(view(x,v1),view(S.nzval,idx),x[col])
        #end
        # idx2 = view(idx,2:length(idx))
        # row = view(S.rowval,idx2)
        # x[row] = view(x,row) .- S.nzval[idx2] .* x[col]
    end
    return x
end

function frw_sb!(x, S, b)
    x .= 0;
    @fastmath @inbounds for col = 1:S.n
        idx = S.colptr[col]+1 : S.colptr[col + 1] - 1
        #println(idx," ",S.rowval[idx])
        x[col] = (b[col] + x[col])/S.nzval[S.colptr[col]]
        #v1 = view(S.rowval,idx)
        @simd for v in idx
             x[S.rowval[v]] -=  S.nzval[v] * x[col]
        end
    end
    return x
end



x[1] = b[1]/S.nzval[1]
x[3] = x[3] - S.nzval[2] * x[1]
x[5] = x[5] - S.nzval[3] * x[1]

x[2] = b[2]/S.nzval[4]
x[5] = x[5] - S.nzval[5]*x[2]

x[3] = (b[3] + x[3])/S.nzval[6]

x[4] = b[4]/S.nzval[7]

x[5] = (b[5] + x[5])/S.nzval[8]
x[6] = x[6] - S.nzval[9]*x[5]

x[6] = (b[6] + x[6])/S.nzval[9]

#--------------------------#
x[6] = b[6]/UU.nzval[9]
x[5] = x[5] - UU.nzval[8] * x[6]

x[5] = (b[5] + x[5]) /UU.nzval[8]
x[2] = x[2] - UU.nzval[7] * x[5]
x[1] = x[1] - UU.nzval[6] * x[5]
