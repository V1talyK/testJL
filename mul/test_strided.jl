using Strided, BenchmarkTools, LinearAlgebra
A = rand(100_000);
B = rand(100_000);
C = rand(100_000);
r = rand(1:100_000,length(A))
r = collect(1:100000)
Bv = view(B,r);
Bv = StridedView(B,r);
@btime copyto!($Br,$Bv);
r1 = CartesianIndices(1:100000)
r2 = CartesianIndices(r)

@btime copyto!($Br,$r1,$B,$r2);
@btime BLAS.blascopy!(100000,$Br,1,$B,1);
@btime $B.=$Bv;

@btime $C.=$A.*$B;
@btime $C.=$A.*$Bv;


@btime @strided $C.=$A.*$B;
@btime @strided $C.=$A.*$Br;

@btime @strided @inbounds permute!($B, $r);

r = collect(1:40)
r = Tuple(r)
