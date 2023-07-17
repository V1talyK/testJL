function get_autoreg(yy)
    ro = zeros(size(yy,1))
    for iw = 1:size(yy,1)
        mu = mean(yy[iw,:])
        s2 = 1/(size(yy,2)-1)*sum(abs2, yy[iw,:] .- mu)
        gam = 1/(size(yy,2)-1)*sum((yy[iw,1:end-1] .- mu).*(yy[iw,2:end] .- mu))
        gam0 = 1/(size(yy,2))*sum((yy[iw,1:end] .- mu).*(yy[iw,1:end] .- mu))

        ro[iw] = gam/gam0
    end
    return ro
end

xx = ones(5,10)
yy = similar(xx)
yy[:,1] .= xx[:,1]
r0 = rand(5)
for t = 2:size(xx,2)
    yy[:,t] .= yy[:,t-1].*r0
end
ro = get_autoreg(yy)
println.(r0, " ", ro)

yy_ar = zeros(33)
pp0 = P0[1]
for t = 1:33
    yy_ar[t] = ro[1]*pp0 + qw[1,t]
    pp0 = yy_ar[t]
end

plt = lineplot(ppl2[1,:])
    lineplot!(plt, yy_ar)
    println(plt)

plt = lineplot(yy_ar)
    println(plt)

test
