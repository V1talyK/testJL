using LinearAlgebra
ACL = cholesky(mA);

@time x0 = ACL\b

@btime L = ACL.L;
@btime LL = LowerTriangular(sparse(L))
Lt = copy(LL');

p = ACL.p;
bp = b[p];
LL[:,p]
@time y .= LL\bp;
@time x .= Lt\y;

sum(abs,x.-x0[p])


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
