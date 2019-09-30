OUT = [joinpath(rdata,"4_mesh_H.tsv"),joinpath(rdata,"4_geom_H.tsv"),joinpath(rdata,"4_wellCon_H.tsv")]

function make_geom(OUT[2])
    
end

function make_mesh(OUT[1])
id = collect(1:length(xy))
xy = collect(Iterators.product(10:20:1000, 10:20:1000))
xy =  convert(Array{Tuple{Float64,Float64},2},xy)
bnd = zeros(Int32,size(xy));
id50 = reshape(id,50,50)
bnd[1,1:end-1] .= id50[1,2:end];
bnd[1:end-1,end] .= id50[2:end,end];
bnd[end,2:end] .= id50[end,1:end-1];
bnd[2:end,1] .= id50[1:end-1,1];

X = map(x->x[1],xy[:]);
Y = map(x->x[2],xy[:]);
bnd = bnd[:]

vxB = Matrix{String}(undef,size(xy))
vyB = Matrix{String}(undef,size(xy))
vxB = map(x->"$(x[1]-10), $(x[1]+10)",xy)
vyB = map(x->"$(x[2]-10), $(x[2]+10)",xy)

vxB[1] = "0, $(vxB[1])"
vxB[50] = "$(vxB[50]), 1000"
vxB[2451] = "0, $(vxB[2451])"
vxB[2500] = "$(vxB[2500]), 1000"

vyB[1] = "0, $(vyB[1])"
vyB[50] = "0, $(vyB[50])"
vyB[2451] = "$(vyB[2451]), 1000"
vyB[2500] = "$(vyB[2500]), 1000"

for i=1:2500
    vxB[i] = "\"$(vxB[i])\""
    vyB[i] = "\"$(vyB[i])\""
end

vxB[2:end,2:end].="\\N"
vyB[2:end,2:end].="\\N"

ioW1 = Base.open(OUT,"w");
writeToPipe(ioW1, id, X, Y, bnd, vxB[:], vyB[:])
close(ioW1)
end
