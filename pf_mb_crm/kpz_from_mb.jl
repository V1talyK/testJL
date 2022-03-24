using XLSX, Dates, SparseArrays, Statistics, UnicodePlots

rfile1 = "/home/lik/Загрузки/Таблица Адаптации из CH.xlsx"
xf = XLSX.readxlsx(rfile1)
sm = XLSX.sheetnames(xf)
DD = xf[sm[1]]
nr = size(DD.dimension)[1]

gi1 = Vector{Int64}(undef,nr-1)
dd1 = Vector{Date}(undef,nr-1)
kpz1 = Vector{Float64}(undef,nr-1)
qz1 = Vector{Float64}(undef,nr-1)

XLSX.openxlsx(rfile1, enable_cache=false) do f
           sheet = f[sm[1]]
           k=0
           k1=0
           for r in XLSX.eachrow(sheet)
              k1+=1
              if k1>2
                  k+=1
              # r is a `SheetRow`, values are read using column references
              rn = XLSX.row_number(r) # `SheetRow` row number
              dd1[k] = Date(r[4])
              kpz1[k] = r[10]
              qz1[k] = r[2]
              gi1[k] = parse(Int64,r[3])
            end
           end
      end

ia = indexin(dd1[2:end-1],vd)
ib = ia[.!isnothing.(ia)]
KPZ = zeros(nt,3)
QZ = zeros(nt,3)
ci = CartesianIndex.(ib,gi1[2:end-1])
KPZ[ci].=kpz1[2:end-1]
QZ[ci].=qz1[2:end-1]

injf = wdd.q.<0
qz = copy(wdd.q)
qz[.!injf].=0

x = []
for i=1:3
    gri1 = indexin(wg[i],wdd.wi)
    gri1 = gri1[.!isnothing.(gri1)]
    push!(x,-KPZ[:,i].*QZ[:,i]/30.5./sum(qz[:,gri1],dims=2)[:])
end

hK1 = deepcopy(hK)
for (k,v) in enumerate(hK1)
    #v["ds"].==vd
    println(k)
    v["kpz"] = isa(v["kpz"],String) ? parse(Float64,v["kpz"]) : v["kpz"]
    v["kpz"] = v["kpz"]*x[in.(v["wi"],wg)][1][Date(v["ds"]).==vd][1]
    v["kpz"] = clamp(v["kpz"],0,1)
end

json_set(r2c,"bondaryK",hK1,cok)
