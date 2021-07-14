well_prod_name = collect("PROD$i" for i = 1:32)
well_inj_name = collect("INJE$i" for i = 1:4)
well_name = vcat(well_prod_name, well_inj_name)

vd = Date(2021,01,01):Dates.Month(1):Date(2030,12,01);

BHPT = zeros(length(well_name),length(vd));

df = DateFormat("dd.mm.yyyy")
D = Vector(undef,0)
for (k,v) in enumerate(Iterators.product(well_name, vd))
    tmp = "$(v[2])"
    push!(D,(v[1], "$(tmp[9:10]).$(tmp[6:7]).$(tmp[1:4])", "BHPT", BHPT[k]))
end

writeToFile("operating_mode_5",[getindex.(D,1),getindex.(D,2),getindex.(D,3),getindex.(D,4)])


#operating_mode_1 - база
#operating_mode_2 - снижение добычи (повышение забойки на доб.)
#operating_mode_3 - снижение закачки (понижение забойки на наг.)
#operating_mode_4 - отключение половины доб. скв.
#operating_mode_5 - отключение наг. скв.
