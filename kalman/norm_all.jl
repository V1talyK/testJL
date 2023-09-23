using UnicodePlots
using Statistics

foo(x) = x^2+x-1
x0 = rand(5:0.1:10, 10000)
x = foo.(x0)

histogram(x, nbins = 11)|>println
lineplot(sort(x))|>println

y = (x.-minimum(x))/(maximum(x) - minimum(x))
histogram(y, nbins = 11)|>println
lineplot(sort(y))|>println

jj = log10(maximum(x) - minimum(x))
y1 = (x.-minimum(x))/10^jj
histogram(y1, nbins = 11)|>println
lineplot(sort(y1))|>println

Q1 = quantile(x,0.25)
Q3 = quantile(x,0.75)
IQR = Q3-Q1

L = Q1 - 1.5*IQR
R = Q3 + 1.5*IQR
y2 = (x .- L)/(R - L)

histogram(y2, nbins = 11)|>println
y3 = copy(y)
for i = 1:1
    a = 3
    y3 = (y3.^a .- 1)./a
end
histogram(y3, nbins = 11)|>println
