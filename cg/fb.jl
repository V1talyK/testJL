function make_CL_in_julia(ACL, nth = 1)
    LL = sparse(ACL.L)
    UU = copy(LL')
    x_temp = zeros(LL.n, nth)
    return (L = LL, U = UU, p = ACL.p, x_temp)
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


function forward_substit!(x, S, b)
    x .= 0;
    @fastmath @inbounds for col = 1:S.n
        idx = S.colptr[col]+1 : S.colptr[col + 1] - 1
        #println(idx," ",S.rowval[idx])
        x[col] = (b[col] + x[col])/S.nzval[S.colptr[col]]
        #v1 = view(S.rowval,idx)
        for v in idx
             x[S.rowval[v]] -=  S.nzval[v] * x[col]
        end
    end
end

function backward_substit!(x, UU, b)
    x .= 0;
    @fastmath @inbounds for col = UU.n:-1:1
        idx = UU.colptr[col + 1]-2 :-1 : UU.colptr[col]
        #println(idx," ",S.rowval[idx])
        x[col] = (b[col] + x[col])/UU.nzval[UU.colptr[col+1]-1]

        for v in idx
             x[UU.rowval[v]] -=  UU.nzval[v] * x[col]
        end
    end
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

x[4] = (b[4] + x[4])/UU.nzval[5]

x[3] = (b[3] + x[3])/UU.nzval[4]
x[1] = x[1] - UU.nzval[3]*x[3]

x[2] = (b[2] + x[2])/UU.nzval[2]

x[1] = (b[1] + x[1])/UU.nzval[1]
