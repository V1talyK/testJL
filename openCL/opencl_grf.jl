using Graphs
g = Graphs.SimpleGraph(3)
add_edge!(g, 1, 2)
add_edge!(g, 3, 2)
add_edge!(g, 3, 1)
plot(g)

using Plots, GraphRecipes
g = wheel_graph(10)
graphplot(g, curves=false)


g = Graphs.SimpleGraph(length(kn[1])+length(kn[2]))

using Clustering
D = rand(100, 100);
D += D'; # symmetric distance matrix (optional)
result = hclust(D, linkage=:single)
cutree(hclu::Hclust; [k], [h])
