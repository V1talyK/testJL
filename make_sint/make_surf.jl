using UnicodePlots
nx = 51
ny = 51

kp = ones(nx,ny);
xy = [[i,j] for i in 1:nx, j in 1:ny][:] # Сетка

xy_kp = collect(zip(rand(1:51,5),rand(1:51,5)))
xy_kp =[(44, 15)
        (28, 26)
        (11, 17)
        (13, 30)
        (43, 43)]
kp0 = [100,250,50,400,100]

mp0 = [0.100,0.250,0.050,0.200,0.100]

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
w_kp = fi0\(kp0.-100)
x = fi0*w_kp .+100

fi = exp.(-(er.*hcat(r...)).^2)
kp = fi'*w_kp .+100

w_mp = fi0\(mp0.-0.1)
x = fi0*w_mp .+0.100
mp = fi'*w_mp .+0.100



plt = heatmap(reshape(kp,nx,ny),colormap=:inferno, height = 51, width = 51);   display(plt);
plt = heatmap(reshape(mp,nx,ny),colormap=:inferno, height = 51, width = 51);   display(plt);

D = reshape(mp,nx,ny)
D = collect(D[i,:] for i=1:nx)
writeToFile("PoroM1",D)

function writeToFile(name,D)
    r=dirname(Base.source_path())
    rfile = joinpath(r,"$name.tsv"); #Синтетика 2
    OUT = joinpath(rfile)
    ioW = Base.open(OUT,"w");
    writeToPipe(ioW, D...)
    close(ioW)
end
