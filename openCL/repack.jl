knLn
fl = repack(knLn, clL, rwL)
ia = findall(fl)

hj = Vector(undef,1)
hj[1] = []
j = 1
upf = false
last_pool = copy(knLn[1])
for (k1, v1) in enumerate(knLn)
    if fl[k1]
        if !upf
            j+=1
            push!(hj,[])
        end
    end
    for (k2,v2) in enumerate(v1)

        push!(hj[j],v2)
        if length(hj[j])==128
             j+=1
             push!(hj,[])
             upf = true
         else
             upf = false
         end
    end
    last_pool = copy(knLn[k1])
end


knLn = hj
vcat(hj...) == vcat(knLn...)

plt = lineplot(length.(hj))
    lineplot!(plt,length.(knLn))
    println(plt)


hj = hj[length.(hj).>0]
findfirst(length.(hj).<128)
hj[97]
hj[98]


findall(length.(hj).=0)

knLn[35]
hj[35]
fl[35:36]





using GraphRecipes, Plots
g = L[1:20,1:20]

graphplot(g, names=1:3, curvature_scalar=0.1)



hj = Vector(undef,6000)
hj = [[] for i in 1:6000]
j=1
last_pool = []
ord = vcat(knLn...)
temp = copy(ord)
temp1 = []
for i=1:30000
    kl = length(temp)
    if kl==0
        break
    end
    v1 = popfirst!(temp)
    c1 = clL[v1]:clL[v1+1]-2
    fl = length(intersect(last_pool,rwL[c1]))==0
    fl1 = issubset(rwL[c1],vcat(hj[1:j-1]...))
    if fl & fl1
        if j==112
            println("$i 1")
            if length(hj[j])==48
                println("hj: ", hj[j])
                println(v1)
                println("48 ", fl)
                println(last_pool)
            end
            if length(hj[j])==1
                println("hj: ", hj[j])
                println(v1)
                println("1 ", fl)
                println(last_pool)
            end
        end
        if length(hj[j])<128
            push!(hj[j],v1)
        else
            j+=1
            push!(hj[j],v1)
        end
    else
        ind = ifelse(length(temp)>=128,128,length(temp)+1)
        if ind>0
            insert!(temp,ind,v1)
        else
            temp = [v1]
        end
    end

    if kl == length(temp)
        j+=1
        if j==112
            println("2")
        end
    end
    last_pool = hj[j]
    if j==112
        println("$i ss")
        if length(hj[j])==49
            println(hj[j])
        end
    end
end
hj = hj[length.(hj).>0]
test_kn(hj)
