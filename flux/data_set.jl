cok=Dict("cma"=>"9c32c492-8265-46f0-a50d-0969d471d303","cdb"=>"","csql"=>"30","ch"=>"ch30"); #Dev 66
pl = "SELECT
iter,
wi,
sm,
sy,
em,
ey,
array(countEqual(st, 0), countEqual(st, 1),countEqual(st, 2), countEqual(st, 3)) as st,
gw,
ws,
kp,
pw,
qw,
ppl,
Jx,
lam,
aqm,
pvt,
nc,
nw,
nt
FROM default._ccord_data_60e400887ddd740040cf4645_ch_scheme_60e40b687ddd740040cf4646 FORMAT JSONCompact"

@time data = chPost(pl,cok)["data"];


x_train = getindex.(data,9)
y_train = getindex.(data,14)

opt = Momentum(0.01, 0.9)

model = Chain(Dense(1377, 500),Dense(500, 8, σ))
model(x_train[1])
loss(x, y) = Flux.Losses.mse(model(x), y)
loss(x_train[1], y_train[1])

data = Vector(undef,450)
for (k,v) in enumerate(zip(x_train, y_train))
    data[k] = Float32.(v[1]), Float32.(v[2])
end

parameters = Flux.params(model)
for i = 1:3
  train!(loss, parameters, data[1:400], opt)
  #println(loss(data[i][1],data[i][2]))
end


mean(loss.(getindex.(data,1),getindex.(data,2)))
loss(data[1][1],data[1][2])
