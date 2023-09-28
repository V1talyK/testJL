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

xx = hcat(y[1:end-1],x[2:end],ones(99))

AA = xx'*xx
bb = xx'*Δy

aa = AA\bb

ac = similar(aa)
ac[1] = aa[1]+1
ac[2] = aa[2]
ac[3] = aa[3]

ax = [aa[2]+1, aa[3], aa[4]]
yx = zeros(100); yx[1] = 1
Δyy = zeros(99)
for i=1:99
    yx[i+1] = ax[1]*yx[i] + ax[2]*x[i+1] + ax[3]

    Δyy[i] = aa[3]*Δx + (1-aa[1])*(y[i-1] - aa[3]/(1-aa[1]) - aa[2]/(1-aa[1])*x[i-1])
end

sum(abs2,y.-yx)

ym = y[1:2:end]
Δym = diff(ym)
xx1 = hcat(ym[1:end-1],x[3:2:end],x[2:2:end-2],ones(49))
AAm = xx1'*xx1
bbm = xx1'*Δym

aam = AAm\bbm
a = sqrt(1+aam[1])
b = aam[2]
b = aam[3]/a
c = aam[4]/(a+1)
