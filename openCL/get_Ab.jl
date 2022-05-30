using JLD, SparseArrays, LinearAlgebra
rsrc1=dirname(dirname(Base.source_path()));
d = JLD.load(joinpath(rsrc1,"myfile.jld"));
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

ACL = cholesky(mA)
L = sparse(ACL.L)
dropzeros!(L)

x0 = L\b
function foo1111(x0)
    x0 .= L\b
end
@btime foo1111($x0)

U = copy(L')
