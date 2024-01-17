aa0 = rand(10)
x = rand(10,100)
y = (aa0'*x)[:] .+ 2

y1, aa1 = lrc(x,y')
mse(vcat(aa0, 2),aa1)

w = ones(11)
w[11] = 0
y2, aa2 = lrcT(x,y'; w = w)
a0 = 20*ones(11)

w3 = 0*ones(11)
w3[11] = 1000
y3, aa3 = lrcGT(x, y'; w = w3, a0 = a0)
