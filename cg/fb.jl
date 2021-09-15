


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
forward_substitution!(LL, bp)
y2 = copy(y)
y1 = copy(y)

@btime x1 = forward_substitution!(UU, y1)
@btime x2 .= forward_substitution1!(UU, y2)
@profiler forward_substitution1!(UU, y2)

sum(abs,x[ACL.p].-x1)

@btime ldiv!($y,$LL,$bp)
@btime ldiv!($x1,$UU,$y)
