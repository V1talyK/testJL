b = rand(20)
A = rand(10,20)
x0 = A'\b
F = svd(A')

M = F.U*diagm(0=>F.S)*F.Vt

M9 = U[:,1:9]*diagm(0=>S[1:9])#*V[:,1:9]'

sum(M.-A')

x9 = M9\b


a0 = rand(10);
x0 = rand(100,10)
y0 = x0*a0
ax0 = (x0'*x0)\(x0'*y0)

F = svd(x0)
x1 = (F.U[:,1:9]*diagm(0=>F.S[1:9])*F.Vt[1:9,:])

ax1 = (x1'*x1)\(x1'*y0)
y1 = x0*ax1

mape(y0,y1)
scatterplot(y0,y1)

rank(x0)
rank(x1)

x0*ax1.-y0
mape(y0,x0*ax1)

xy = collect(Iterators.product(0:100,0:100))
x = getindex.(xy,1)
y = getindex.(xy,2)
z = x/20 .^2 .+ sin.(y)./5 .+ y/10 .^2 .+ x.*(y.-20)/1000
z[40 .<x .<60] .= 5
z[20 .<y .<25] .= -1
    plt = heatmap(z, height=50);
    println(plt)

F = svd(z)
z1 = (F.U[:,1:2]*diagm(0=>F.S[1:2])*F.Vt[1:2,:])

plt = heatmap(z1, height=50);
    println(plt)

mse(z, z1)
