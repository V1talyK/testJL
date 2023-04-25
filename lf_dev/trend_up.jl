using UnicodePlots, StatsBase
x0 = collect(1:0.1:10)
y0 = 2 .*x0 .+ 1



y1 = 2 .*x0 .+ 5
y1 = y0 + 1*sin.(x0*5)
lineplot(y1)|>println
y1[5] = 7

LF = sum(abs2,mean(abs, y0 .- y1) .- abs.(y0 .- y1))
LFf(y1) = sum(abs2,mean(y0 .- y1) .- (y0 .- y1))

fg = Vector(undef, 10)
for i=1:10
    y1[i] = y1[i] - 5
    fg[i] = LFf(y1)
end
