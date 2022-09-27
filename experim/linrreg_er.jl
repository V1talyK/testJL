using UnicodePlots, StatsBase, HypothesisTests
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

N = length(x)
Sy2 = 1/(N-2)*(sum(y.^2) - a1*sum(x.*y)-b1*sum(y))
Sa2 = Sy2/(sum(x.^2) - N*x1^2)
Sb2 = Sy2*sum(x.^2)/(N*sum(x.^2) - N*x1^2)

tp = 2.13
Δa = tp*sqrt(Sa2)
Δb = tp*sqrt(Sb2)

"$a1+$Δa"
"$b1+$Δb"

foo_dy(x) = sqrt(x^2*Δa*Δa + 1 * Δb*Δb)
Δy = foo_dy.(x)
plt = histogram(Δy)
    println(plt)

plt = scatterplot(x,Δy)
    println(plt)
