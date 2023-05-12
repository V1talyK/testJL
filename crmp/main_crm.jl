using UnicodePlots
qp0 = rand(2:5,5,10)
qp = rand(2.0:5.0,5,10)
qi = hcat(fill(rand(2:5,4),10)...)
qi[1,6:10] .= 5
pw = 10*ones(5,11)
dpw = diff(pw, dims=2)
bet = (3.7e-4 + 7.4e-3)/2;
tay = 250/3*250/3*1*0.14*bet/0.5#100
fij = ones(4,5)./5
ei = 0.
Jj = 1
dt = 30.5
nt = 10
tk = range(0; length = nt, step = 30.5)

for j=1:5
    for t = 1:nt
        ss = 0.0
        for k=1:t
            ss+=(ei + sum(fij[:,j].*qi[:,t]) - Jj*tay*dpw[j,k]/dt)*exp((tk[k] - tk[t])/tay)*(1-exp(-dt/tay))
        end
        qp[j,t] = qp0[j,1]*exp((0-tk[t])/tay) + ss
    end
end

lineplot(qp[2,:]) |> println
lineplot(qi[1,:]) |> println
