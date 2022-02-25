using LinearAlgebra, BenchmarkTools
LinearAlgebra.BLAS.set_num_threads(1)

r, c, v = findnz(A)
dP = x[r].-x[c]
@btime acP = accumarray(r,dP)
xr = view(x,r)
xc = view(x,c)
@btime begin
    @btime @inbounds dP .= xr.-xc
    @btime accumarray!(acP, r, dP)
end


A1 = sparse(r,c,1)
A1 = A1 .- sparse(1:size(A,1),1:size(A,1),accumarray(r,ones(length(r))))

@btime dP2 .= A1*x
@btime mul!(dP2, A1, x)


sum(abs,acP.+dP2)
