function unpackSh(sh)

    nrow = size(sh)[1]
    wi = zeros(Int64,nrow)
    wn = Vector{String}(undef,nrow)
    d = Vector{Date}(undef,nrow)
    st = zeros(Int64,nrow)
    ufl = zeros(Int64,nrow)
    pw = zeros(nrow)
    qw = zeros(nrow)
    ws = zeros(nrow)
    opra = zeros(nrow)
    fl = falses(nrow)

    for i=1:nrow
        #println(i)
        fl[i] = !ismissing(sh[i,1]) | !ismissing(sh[i,2])
        wi[i] = isa(sh[i,1],String) ? tryparse(Int64,sh[i,1]) : ifelse(ismissing(sh[i,1]),0,sh[i,1])
        wn[i] = ismissing(sh[i,2]) ? "" : sh[i,2]
        d[i] = ismissing(sh[i,3]) ? Date("2020-01-01") : date_iso(sh[i,3])
        st[i] = ismissing(sh[i,4]) ? 0 : ifelse(sh[i,4]=="Раб.",1,0)
        ufl[i] = ismissing(sh[i,5]) ? 0 : ifelse(sh[i,5]=="Деб.",1,0)
        pw[i] = ismissing(sh[i,6]) ? NaN : sh[i,6]
        qw[i] = ismissing(sh[i,7]) ? NaN : isa(sh[i,7],String) ? tryparse(Float64,sh[i,7]) : sh[i,7]
        ws[i] = ismissing(sh[i,8]) ? NaN : isa(sh[i,8],String) ? tryparse(Float64,sh[i,8]) : sh[i,8]
        opra[i] = ismissing(sh[i,9]) ? NaN : isa(sh[i,9],String) ? tryparse(Float64,sh[i,9]) : sh[i,9]
    end

    wi = wi[fl];
    wn = wn[fl];
    d = d[fl];
    st = st[fl];
    ufl = ufl[fl];
    pw = pw[fl];
    qw = qw[fl];
    ws = ws[fl];
    opra = opra[fl];
    return wi, wn, d, st, ufl, pw, qw, ws, opra
end

function chPost(pl,cok)
   prx=getProxy();
   return chPost(pl,cok,prx)
end

function chPost(pl,cok,prx)
    ch = cok["ch"];
    # res=Requests.post("http://$prx/$ch/", data = pl)
    # res = readstring(res);
    #pl = JSON.json(pl)
    #println(pl)
    headers = Dict("Cookie" => "ma_session = $(cok["cma"]); ma_db = $(cok["cdb"])")
    if haskey(cok,"Authorization")
        headers["Authorization"]=cok["Authorization"]
    end
    #println("http://$prx/$ch/")
    res = HTTP.request("POST", "http://$prx/$ch/", headers, pl);
    if res.status!=200
        #if res[1:4]=="Code"
            println(res[1:100])
        #end;
    end;
    # println("----------")
    # println(String(res.body)=="")
    # println("----------")
    Sres = String(res.body);
    if Sres!=""
        return JSON.parse(Sres)
    else
        return Sres
    end
end

function setUTlb2CH(chname, cok, rcd::Vector{Tuple})
    #Сохраняем сохраняем универсальный json в ClickHouse по имени из вне
    name, dbname = buildName(chname);
  return funSetUTlb2CH(rcd,dbname,name,cok)
end

function funSetUTlb2CH(rcd::Vector{Tuple},dbname,name,cok)
    rcd = replace.(string.(rcd),"\""=>"")
    rcd = replace.(rcd,"Int32"=>"")
    rcd = replace.(rcd,"Int8"=>"")
    rcd = replace.(rcd,"Float32"=>"")
    rcd = replace.(rcd,"Any"=>"")
    rcd = replace.(rcd,"Vector{}"=>"")
    pl = string("INSERT INTO $dbname.$name  values ",join(rcd," "));
    #println(pl)
    res = chPost(pl,cok);
    #pl = string("INSERT INTO $dbname.$name  values ");
    #res = chStreamPost(vcat(pl,p2),cok)
   return res
end

function getWells(ofid::Int64,cok::Dict)
    res = sqlGet(cok,"well?select=name,id&oil_field_id=eq.$ofid")
    wi = Base.get.(res,"id",0)
    wn = Base.get.(res,"name","")
    return wi, wn
end

function sqlGet(cok,s3ng)
    prx=getProxy();
    sql = haskey(cok,"sql") ? cok["sql"] : "sql$(cok["csql"])";
    #res = Requests.get("http://$proxy/$sql/$s3ng")
    headers = Dict()
    if haskey(cok,"Authorization")
        headers["Authorization"]=cok["Authorization"]
    end

    res = HTTP.request("GET", "http://$prx/$sql/$s3ng",headers)
    res = JSON.parse(String(res.body));
    return res
end

function getProxy()
  proxy = isdefined(Main, :proxy) ? Main.proxy : "proxy"
  return proxy
end

function get_object(ref,cok="123")
    data = JSON.json(Dict("list"=>[ref]))
    res = mongoPost("get_objects", data, cok)
  return res["data"]["result"][1]
end

function makeRef(ref::String,id::String)
    return r2c=Dict("\$ref"=>ref,"\$id"=> Dict("\$oid"=> id))
end
function makeRef(s3ng::String)
    s3ng = split(s3ng,":")
    return makeRef(String(s3ng[1]),String(s3ng[2]))
end
function buildName(dbname::String, v::Vector)
    v = map(x->string("_",x["\$ref"],"_",x["\$id"]["\$oid"]),v)
    name=string(dbname,".",string(v...))
end

function mongoPost(fun, data, cok)
  prx=getProxy();
  return mongoPost(fun, data, cok, prx)
end

function mongoPost(fun, data, cok, prx)
  ma = haskey(cok,"ma") ? cok["ma"] : "ma";
  #cks = Dict("ma_session" => cok["cma"], "ma_db" => cok["cdb"])

  headers = Dict("Cookie" => "ma_session = $(cok["cma"]); ma_db = $(cok["cdb"])")
  if haskey(cok,"Authorization")
      headers["Authorization"]=cok["Authorization"]
  end
  res = HTTP.request("POST", "http://$prx/$ma/$fun", headers, data);
  #res = Requests.post("http://$prx/$ma/$fun", data = data, cookies = cks)
  #res=readstring(res);
  res = JSON.parse(String(res.body));
  #res=JSON.parse(res);
  return res
end

function date_iso(diso::String)
    #Ускоренное преобразование панорамных дат формата YYYY-MM-DD
    year = parse(Int, diso[1:4])
    month = parse(Int, diso[6:7])
    day = length(diso)>=10 ? parse(Int, diso[9:10]) : 1;
    return Date(year, month, day)
end
