using  UnicodePlots, StatsBase, Term
include(joinpath(Base.source_path(),"../../sim/sim9.jl"))


nw, nt  = 9, 100;
qw = rand(nw,nt)
qw[1,:] .= 1;
qw[2,:] .= 1;
qw[3,:] .= 1;
qw[4,:] .= 1;
qw[5,:] .= -8; qw[5,25:50] .= -6;
qw[6,:] .= 1;
qw[7,:] .= 1;
qw[8,:] .= 1;
qw[9,:] .= 1;

Tt = rand(nw)

Pa = 10;
P0 = Pa*ones(nw)
bet = 1e-4;

AA, bb, eVp = make9p(Tt, Pa, nw, bet)

PM = sim(qw, nt, AA, bb, P0, eVp)
plt = Vector(undef, nw)
for i=1:nw
  plt[i] = lineplot(PM[i,:])
end

grid(panel.(plt[1:3]); layout=(1, nothing)) |> print
grid(panel.(plt[4:6]); layout=(1, nothing)) |> print
grid(panel.(plt[7:9]); layout=(1, nothing)) |> print
