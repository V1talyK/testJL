A = rand(21,21)
A = (A.+A')./2
A[1:size(A,1)+1:end].=0

vb = Vector(undef, size(A,1))
AA = A.>=0.75
sp = findall(any(AA, dims = 2)[:])

for i in sp
    vb[i] = []
    ia = findall(AA[i,:])
    for j = 1:length(ia)
        vp = [i]
        flag = true
        fg = ia[j]
        while flag
            vp = vcat(vp, fg)
            ib = findall(prod(AA[vp,:], dims=1)[:])
            flag = length(ib)>0
            if flag
                fg = ib[1]
            else
                fg = []
            end
        end
        push!(vb[i],vp)
    end
end

for (k, v) in enumerate(vb)
    vb[k] = unique(Set.(v))
end

vb = unique(vcat(vb...))
