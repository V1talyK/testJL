using UnicodePlots, StatsBase
rsrc = dirname(Base.source_path())
infile = string(rsrc,"/9w.csv")
DT = Vector(undef,0)
open(infile, "r") do f
      linecounter = 0
      timetaken = @elapsed for l in eachline(f)
         linecounter += 1
         if linecounter>1
            push!(DT,split(l,","))
         end
      end
      (timetaken, linecounter)
 end

wi = tryparse.(Int64,getindex.(DT,1))
ti = tryparse.(Float64,strip.(getindex.(DT,3),'"')); ti = round.(Int64, ti)
ql = tryparse.(Float64,strip.(getindex.(DT,4),'"'))
qo = tryparse.(Float64,strip.(getindex.(DT,5),'"'))
qw = tryparse.(Float64,strip.(getindex.(DT,6),'"'))
qi = tryparse.(Float64,strip.(getindex.(DT,7),'"'))
pw = tryparse.(Float64,strip.(getindex.(DT,8),'"'))
ppl = tryparse.(Float64,strip.(getindex.(DT,9),'"'))

uwi = sort(unique(wi)); nw = length(uwi);
uti = sort(unique(ti)); nt = length(uti);

r = indexin(wi,uwi)
c = indexin(ti,uti)

ql1 = zeros(nw, nt)
qo1 = zeros(nw, nt)
qw1 = zeros(nw, nt)
qi1 = zeros(nw, nt)
pw1 = zeros(nw, nt)
ppl1 = zeros(nw, nt)

for (k,v) in enumerate(zip(r,c))
   ql1[v[1],v[2]] = ql[k]
   qo1[v[1],v[2]] = qo[k]
   qw1[v[1],v[2]] = qw[k]
   qi1[v[1],v[2]] = qi[k]
   pw1[v[1],v[2]] = pw[k]/10
   ppl1[v[1],v[2]] = ppl[k]/10
end

lineplot(ql1[5,:]) |> println
lineplot(pw1[5,:]) |> println

ql2 = zeros(nw, round(Int64,ceil(nt/30)))
qo2 = zeros(nw, round(Int64,ceil(nt/30)))
qw2 = zeros(nw, round(Int64,ceil(nt/30)))
qi2 = zeros(nw, round(Int64,ceil(nt/30)))
pw2 = zeros(nw, round(Int64,ceil(nt/30)))
ppl2 = zeros(nw, round(Int64,ceil(nt/30)))

for (k,v) in enumerate(Iterators.partition(uti,30))
      ql2[:,k] .= mean(ql1[:,v],dims=2)
      qo2[:,k] .= mean(qo1[:,v],dims=2)
      qw2[:,k] .= mean(qw1[:,v],dims=2)
      qi2[:,k] .= mean(qi1[:,v],dims=2)
      pw2[:,k] .= mean(pw1[:,v],dims=2)
      ppl2[:,k] .= mean(ppl1[:,v],dims=2)
end

lineplot(ql2[5,:]) |> println
lineplot(pw2[5,:]) |> println

lineplot(qi2[3,:]) |> println

iw = 8
   lineplot(ql2[iw,:]./(ppl2[iw,:].-pw2[iw,:])) |> println
