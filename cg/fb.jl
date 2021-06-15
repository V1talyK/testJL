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

function forward_substitution!(L, y::AbstractVector)
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

@btime $x .= $CL\$b

LL = LowerTriangular(sparse(CL.L))
bp = copy(b[CL.p])
UU = copy(transpose(LL))


@btime y = backward_substitution!(LL, bp)
y2 = copy(y)
y1 = copy(y)

@btime x1 = forward_substitution!(UU, y1)
@btime $x2 = forward_substitution1!($UU, $y2)
@profiler forward_substitution1!(UU, y2)

sum(abs,x[CL.p].-x1)

@btime ldiv!($y,$LL,$bp)
@btime ldiv!($x1,$UU,$y)
