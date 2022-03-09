using UnicodePlots, Statistics, Dates
nx = 51
ny = 51

kp = ones(nx,ny);
xy = [[i,j] for i in 1:nx, j in 1:ny][:] # Сетка

xy_kp = collect(zip(rand(1:51,10),rand(1:51,10)))
xy_kp =[(44, 15)
        (28, 26)
        (11, 17)
        (13, 30)
        (43, 43)]
# xy_kp =[(51, 32)
        # (6, 21)
        # (6, 29)
        # (7, 49)
        # (45, 13)
        # (39, 37)
        # (47, 14)
        # (26, 22)
        # (12, 29)
        # (15, 20)]

#kp0 = rand(25:200,10);#
kp0 = [100,250,50,400,100]
#mp0 = rand(15:25,10)./100;
mp0 = [0.15,0.250,0.20,0.200,0.15]

r = Vector(undef, length(xy))
for (k, v) in enumerate(xy)
    #r[k] = sqrt(sum((v.-xy_kp').^2))
    r[k]= sqrt.((getindex.(xy_kp,1) .- v[1]).^2 .+ (getindex.(xy_kp,2) .- v[2]).^2)
end

r0 = Vector(undef, length(xy_kp))
for (k, v) in enumerate(xy_kp)
    r0[k]= sqrt.((getindex.(xy_kp,1) .- v[1]).^2 .+ (getindex.(xy_kp,2) .- v[2]).^2)
end

er = 1e-1
fi0 = exp.(-(er.*hcat(r0...)).^2)
w_kp = fi0\(kp0.-median(kp0))
x = fi0*w_kp .+median(kp0)

fi = exp.(-(er.*hcat(r...)).^2)
kp = fi'*w_kp .+median(kp0)

w_mp = fi0\(mp0.-median(mp0))
x = fi0*w_mp .+median(mp0)
mp = fi'*w_mp .+median(mp0)



plt = heatmap(reshape(kp,nx,ny),colormap=:inferno, height = 51, width = 51);   display(plt);
plt = heatmap(reshape(mp,nx,ny),colormap=:inferno, height = 51, width = 51);   display(plt);

D = reshape(mp,nx,ny)
D = collect(D[i,:] for i=1:nx)
writeToFile("export/PoroM2",D)
