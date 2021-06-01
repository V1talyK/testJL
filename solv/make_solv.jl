using LinearAlgebra, BenchmarkTools
LinearAlgebra.BLAS.set_num_threads(1)

ACL = cholesky(mA);

@btime x0 = ACL\b

@btime L = ACL.L;
@btime LL = LowerTriangular(sparse(L))
Lt = copy(LL');

p = ACL.p;
bp = b[p];
LL[:,p]
@time y .= LL\bp;
@time x .= Lt\y;

sum(abs,x.-x0)


@btime ldiv!(y, LL,bp)
@btime ldiv!(x, Lt,y)
@btime x0 = ACL\b

@btime begin
    ldiv!(y, LL,bp)
    ldiv!(x, Lt,y)
end

LinearAlgebra.BLAS.trsv('L', 'N', 'N', LL, bp)

LinearAlgebra.BLAS.set_num_threads(2)

@profiler ldiv!(y, LL,bp)


permute(sparse(L), collect(1:16610), p)


function mkslv(ACL)
    LL = LowerTriangular(sparse(ACL.L))
    Lt = copy(LL');
    p = ACL.p;
    invp = invperm(p)
    x = zeros(length(p))
    y = zeros(length(p))
    function slvr!(_x,b)
        bp = view(b,p);
        ldiv!(y, LL,bp)
        ldiv!(x, Lt,y)
        return _x[p] = x
    end
end
function mkslv(ALU)
    LL = LowerTriangular(sparse(ALU.L))
    UU = UpperTriangular(sparse(ALU.U))
    p = ALU.p;
    invp = invperm(p)
    x = zeros(length(p))
    y = zeros(length(p))
    function slvr!(_x,b)
        bp = view(b,p);
        ldiv!(y, LL,bp)
        ldiv!(x, UU,y)
        return _x[p] = x
    end
end

slvr! = mkslv(ALU)
x = zeros(length(b))
@btime slvr!(x,b)

LL = LowerTriangular(sparse(ALU.L))
UU = LowerTriangular(sparse(ALU.U))
@btime ldiv!(x, UU, b)
ALU = lu(mA)
