using Clustering
using Plots
using SparseArrays
D = rand(10, 10)
D += D'
hc = hclust(D, linkage=:single)
plot(hc)

AA = rand(Float64, 100,100);
AA.+=AA';
AA[1:101:end].+=100

AA[1,end]  = 0
sAA = sparse(AA)
@btime exp($AA);
@btime inv($AA);
@btime $AA*$AA;