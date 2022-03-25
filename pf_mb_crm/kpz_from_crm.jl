using XLSX, Dates, SparseArrays, Statistics, UnicodePlots

rfile1 = "/media/lik/3808EADB08EA96E2/YaDisk/Works/ТатНефть/Ромашкинское/результаты/вр/КПЗ по CRM - Сводный файл.xlsx"
xf = XLSX.readxlsx(rfile1)
sm = XLSX.sheetnames(xf)
DD = xf[sm[1]]
nr = size(DD.dimension)[1]

wi1 = Vector{Int64}(undef,nr-1)
dd1 = Vector{Date}(undef,nr-1)
kpz1 = Vector{Float64}(undef,nr-1)
qz1 = Vector{Float64}(undef,nr-1)

XLSX.openxlsx(rfile1, enable_cache=false) do f
           sheet = f[sm[1]]
           k=0
           k1=0
           for r in XLSX.eachrow(sheet)
              k1+=1
              if k1>1
                  k+=1
              # r is a `SheetRow`, values are read using column references
              rn = XLSX.row_number(r) # `SheetRow` row number
              dd1[k] = Date(r[3])
              kpz1[k] = r[9]
              qz1[k] = r[4]
              wi1[k] = parse(Int64,r[2])
            end
           end
      end

ia = indexin(dd1,vd)
ib = ia[.!isnothing.(ia)]
KPZ = zeros(nw,nt)
QZ = zeros(nw,nt)
ic = indexin(wi1,wdd.wi)[.!isnothing.(ia)]
ci = CartesianIndex.(ic,ib)
KPZ[ci].=kpz1
QZ[ci].=qz1

injf = wdd.q.<0
qz = copy(wdd.q)
qz[.!injf].=0

sum(qz'[QZ.!=0].+QZ[QZ.!=0]/30.5)

x = []
for i=1:3
    gri1 = indexin(wg[i],wdd.wi)
    gri1 = gri1[.!isnothing.(gri1)]
    push!(x,-KPZ[:,i].*QZ[:,i]/30.5./sum(qz[:,gri1],dims=2)[:])
end


hK1 = deepcopy(hK)
for (k,v) in enumerate(hK1)
    #v["ds"].==vd
    #println(k)
    v["kpz"] = isa(v["kpz"],String) ? parse(Float64,v["kpz"]) : v["kpz"]
    t1 = kpz1[Date(v["ds"]).==dd1]
    t2 = wi1[Date(v["ds"]).==dd1]
    t3 = t1[t2.==v["wi"]]
    if length(t3)>0
        println(k," ",v["kpz"]," ",t3[1])
        v["kpz"] = t3[1]
    end
    #
    v["kpz"] = clamp(v["kpz"],0,1)
end

json_set(r2c,"bondaryK",hK1,cok)
