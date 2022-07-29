using UnicodePlots, StatsBase
a = 0.5
b = 2.5
foo(x) = a*x+b

x = rand(0:10,16)
y = foo.(x) .+ rand(-1:0.125:1,length(x))

plt = scatterplot(x,y)
    lineplot!(plt,x,foo.(x))
    println(plt)

x1 = mean(x)
y1 = mean(y)
a1 = sum((x.-x1).*(y.-y1))/sum((x.-x1).^2)
b1 = y1 - a1*x1

lineplot!(plt,x,a1.*x.+b1)
    println(plt)


Sy2 = 1/(length(x)-2)*(sum)
