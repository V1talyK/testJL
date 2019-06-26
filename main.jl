using JLD, SparseArrays
r=dirname(Base.source_path());
d = JLD.load(joinpath(r,"myfile.jld"));
d = JLD.load(joinpath(r,"myfile200k.jld"));
A=d["A"];

c = Vector(undef,length(A.rowval))
for i=1:length(A.colptr)-1
    ind = convert(Array{Int64,1},A.colptr[i]:(A.colptr[i+1]-1))
    c[ind] = i*ones(Int64,length(ind));
end
A = SparseArrays.sparse(c, A.rowval, A.nzval)
#A = convert(SparseArrays.SparseMatrixCSR{Float64},A)
b=d["b"];
mA = -A;


@time x = A\b;
using AlgebraicMultigrid
@time begin
    for i=1:100
        ml = ruge_stuben(A);
        AlgebraicMultigrid.solve(ml, b);
    end
end

@time p = aspreconditioner(ml)
@time c = cg(A, b, Pl = p)
@time solve(ml, b, tol=1e-5);

@time x1 = mA\b;

@time LUf=LinearAlgebra.lu(A);
@time for i=1:1
    x3=LUf\b;
end
@time CL=cholesky(-A);
@time for i=1:1
    x2 = CL\b;
end

L = LU[:L];
Rs = LU[:Rs];
p = LU[:p];
U = LU[:U];
@time y=L\view(Rs.*b,p);
@time x1=U\y;
@time A_ldiv_B!(x,LU,b)
@time x1 = A_ldiv_B!(CL,b)
y = CL[:L]\b
x1 = CL[:L]'\y
x[1]
x1[LU[:p].==1][1]

sum(LU[:L]*LU[:U]-(LU[:Rs] .* A)[LU[:p], LU[:q]])


@time LU[:L]

Y = zeros(Float64,length(b));
@time A_ldiv_B!(Y, LU, b)
