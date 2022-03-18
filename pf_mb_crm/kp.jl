using XLSX, Dates, SparseArrays, Statistics, UnicodePlots

rfile1 = "/media/lik/3808EADB08EA96E2/YaDisk/Works/ТатНефть/Ромашкинское/результаты/Общие результаты_КПЗ.xlsx"
xf = XLSX.readxlsx(rfile1)
sm = XLSX.sheetnames(xf)
DD = xf[sm[1]]
nr = size(DD.dimension)[1]

wi1 = Vector{Int64}(undef,nr-1)
wn1 = Vector{String}(undef,nr-1)
dd1 = Vector{Date}(undef,nr-1)
kpz1 = Vector{Float64}(undef,nr-1)

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
              wn1[k] = isa(r[2],Int64) ? string(r[2]) : r[2]
              dd1[k] = Date(r[3])
              kpz1[k] = r[10]
              wi1[k] = parse(Int64,r[1])
            end
           end
      end


rfile2 = "/media/lik/3808EADB08EA96E2/YaDisk/Works/ТатНефть/Ромашкинское/результаты/Сводный файл по CRM.xlsx"
xf = XLSX.readxlsx(rfile2)
sm = XLSX.sheetnames(xf)
DD = xf[sm[2]]
nr = size(DD.dimension)[1]

wi2 = Vector{Int64}(undef,nr-1)
wn2 = Vector{String}(undef,nr-1)
dd2 = Vector{Date}(undef,nr-1)
kpz2 = Vector{Float64}(undef,nr-1)

XLSX.openxlsx(rfile2, enable_cache=false) do f
         sheet = f[sm[2]]
         k=0
         k1=0
         for r in XLSX.eachrow(sheet)
            k1+=1
            if k1>2
                k+=1
            # r is a `SheetRow`, values are read using column references
            rn = XLSX.row_number(r) # `SheetRow` row number
            wn2[k] = isa(r[1],Int64) ? string(r[1]) : r[1]
            dd2[k] = Date(r[3])
            kpz2[k] = r[9]
            wi2[k] = parse(Int64,r[2])
          end
         end
    end

    uwi1 = unique(wi1)
    uwi2 = unique(wi2)

    uwi = intersect(uwi1,uwi2)

    wi1,kpz1,dd1,wn1 = flt_by_wi(uwi,wi1,kpz1,dd1,wn1)
    wi2,kpz2,dd2,wn2 = flt_by_wi(uwi,wi2,kpz2,dd2,wn2)

    evd = extrema(dd1)
    evd2 = extrema(dd2)
    vd = evd2[1]:Dates.Month(1):evd[2]

    wi1,kpz1,dd1,wn1 = flt_by_dd(vd,wi1,kpz1,dd1,wn1)
    wi2,kpz2,dd2,wn2 = flt_by_dd(vd,wi2,kpz2,dd2,wn2)

    ib1 = indexin(dd1,vd)
    ib2 = indexin(dd2,vd)

    ia1 = indexin(wi1,uwi)
    ia2 = indexin(wi2,uwi)

    KPZ1 = zeros(length(uwi),length(vd))
    KPZ2 = zeros(length(uwi),length(vd))

    ci1 = CartesianIndex.(ia1,ib1)
    ci2 = CartesianIndex.(ia2,ib2)

KPZ1[ci1].=kpz1
KPZ2[ci2].=kpz2

KPZ11 = KPZ1[.&(KPZ1.!=0,KPZ2.!=0)]
KPZ22 = KPZ2[.&(KPZ1.!=0,KPZ2.!=0)]

plt = scatterplot(KPZ11,KPZ22)
    println(plt)

mean(abs2.(KPZ11.-KPZ22))
mean(abs.((KPZ11.-KPZ22./KPZ11)))
cor(KPZ11,KPZ22)


function flt_by_wi(uwi,wi,kpz,dd,wn)
    ia = indexin(wi,uwi)
    ib = .!isnothing.(ia)
    return wi[ib], kpz[ib], dd[ib], wn[ib]
end

function flt_by_dd(vd,wi,kpz,dd,wn)
    ia = indexin(dd,vd)
    ib = .!isnothing.(ia)
    return wi[ib], kpz[ib], dd[ib], wn[ib]
end
