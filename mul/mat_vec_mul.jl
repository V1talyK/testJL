using BenchmarkTools, LoopVectorization

x1 = similar(x)
x2 = similar(x)
x3 = similar(x)
@btime foo($x1, $A, $b)
@btime smv!($x1, $A, $b)
@btime acar!($x1, $r,$c,$v, $b)

foo(x1, A, b)
smv!(x2, A, b)
acar!(x3, r,c,v, b)

sum(abs,x1.-x2)
sum(abs,x1.-x3)

@profiler smv!(x1, A, b)

function foo(x,A,b)
    x .= A*b
end

function smv!(x, A, b)
    @inbounds @fastmath for k = 1:length(x)
        s = zero(eltype(x))
        @simd for kk = A.colptr[k]:A.colptr[k+1]-1
            s+=A.nzval[kk]*b[A.rowval[kk]]#*-1
        end
        x[k] = s
    end
end

r,c,v = findnz(A)

function acar!(x, r, c, v, b)
    x.=0.0
    @turbo for i in eachindex(r)
        x[r[i]]+=v[i]*b[c[i]]
    end
end
