using JLD
r=dirname(Base.source_path());
d = JLD.load(joinpath(r,"myfile.jld"));
A=d["A"];
b=d["b"];


@time x= A\b;

@time LU=lufact(A);
L = LU[:L];
Rs = LU[:Rs];
p = LU[:p];
U = LU[:U];
@time y=L\view(Rs.*b,p);
@time x1=U\y;

x[1]
x1[LU[:p].==1][1]

sum(LU[:L]*LU[:U]-(LU[:Rs] .* A)[LU[:p], LU[:q]])


@time LU[:L]

Y = zeros(Float64,length(b));
@time A_ldiv_B!(Y, LU, b)
