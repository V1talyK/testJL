using UnicodePlots

y = zeros(100)
x = rand(100)
a = [0.5, 2, 4]
y[1] = 1

for i=1:99
    y[i+1] = a[1]*y[i] + a[2]*x[i+1] + a[3]
end

lineplot(y)|>println

Δy = diff(y)
Δx = diff(x)

#Δy = a[2]*Δx - (1-a[1])*(y[i-1] - a[3]/(1-a[1]) - a[2]/(1-a[1])*x[i-1])

xx = hcat(Δx,y[1:end-1],x[1:end-1],ones(99))

AA = xx'*xx
bb = xx'*Δy

aa = AA\bb

ac = similar(aa)
ac[2] = aa[1]
ac[1] = aa[2]+1
ac[3] = aa[4]

ax = [aa[2]+1, aa[3], aa[4]]
yx = zeros(100); yx[1] = 1
Δyy = zeros(99)
for i=1:99
    yx[i+1] = ax[1]*yx[i] + ax[2]*x[i+1] + ax[3]

    Δyy[i] = aa[3]*Δx + (1-aa[1])*(y[i-1] - aa[3]/(1-aa[1]) - aa[2]/(1-aa[1])*x[i-1])
end

sum(abs2,y.-yx)
