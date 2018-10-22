using JLD
r=dirname(Base.source_path());
d = JLD.load(joinpath(r,"myfile.jld"));
A=d["A"];
b=d["b"];
mA = -A;

@time x = A\b;
@time x1 = mA\b;

@time LU=lufact(A);
@time for i=1:50
    x3=LU\b;
end
@time CL=cholfact(-A);
@time for i=1:50
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
