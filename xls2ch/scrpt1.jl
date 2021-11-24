using XLSX, Dates
include(joinpath(dirname(Base.source_path()),"libs.jl"))


r2c=makeRef("ccord_data:61691e5a7ddd74004edb10b4");#Тест
r2c=makeRef("ccord_data:61697b377ddd74004edb1103");#Елабуга_Башкир_Анализ разработки
scheme=makeRef("ch_scheme:6059ccb67ddd740041ed9e07");

obj=get_object(r2c, cok)["items"];
expv = obj["opt_link"]
chname = buildName("default",[expv,scheme]);
ofid = get_object(get_object(obj["dev_object"],cok)["items"]["oilfield"],cok)["items"]["code"]

wi0, wn0 = getWells(ofid,cok)

rxf = XLSX.readxlsx("/home/lik/proto/testJL/xls2ch/exp_v7.xlsx")
shlist = XLSX.sheetnames(rxf)
sh = rxf[XLSX.sheetnames(rxf)[1]]

sh = XLSX.readdata("/home/lik/proto/testJL/xls2ch/exp_v7.xlsx", shlist[1], "A2:I93")
wi = getindex.(sh,1)
wi, wn, d, st, ufl, pw, qw, ws, opra = unpackSh(sh);

ia = indexin(wn, wn0);
wi = wi0[ia]

#genUTlb4CH(chname, cok, clm)

rcdB = Vector{Tuple}(undef, count(fl));
k = 0
for i=1:count(fl)
    rcdB[i] = (wi[i],string("'",d[i],"'"),st[i], ufl[i], pw[i], 0, qw[i], 0, ws[i],opra[i])
end
setUTlb2CH(chname,cok, rcdB)
