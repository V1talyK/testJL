ACR = makeUniqTimeInterv(flagUpACL,3)

makeDBofGrad(flagUpACL)

uACR = unique(ACR[:,3:4],dims=1);
B1 = falses(size(uACR,1));
B2 = zeros(size(uACR,1));
D =  zeros(nt,3);
k = [0]
ui = CartesianIndex.(uACR[:,1],uACR[:,2])


for t=1:nt
    prox2mat!(ACR,t,ui,B1,B2,(x->x))
end


fc = Vector(undef, nt)
fg = Vector(undef, nt)
fp = Vector(undef, nt)
nl = 3

MB = zeros(nc+nw,nw*nl)
flag = falses(nw*nl)
ACCL = Vector(undef,7);
for i=1:7
    ACCL[i] = spzeros(nc+nw,nc+nw);
    updatesp!(ACCL[i],1:nc+nw,1:nc+nw,-i)
end
B = zeros(nc+nw);
HJ = Matrix(undef,nt,nl);
for i=1#:nw
    for t=1:nt-3
        for j=1:nl
            ia = findfirst(findlast(flagUpACL[1:t+j-1]).==findall(flagUpACL))
            HJ[t,j] = CartesianIndex(ia,j);
            if CartesianIndex(ia,j) in HJ[1:t-1,:]

            else
                B.=0;
                if j == 1
                    B[nc+1:nc+nw].=1;
                else
                    B[1:nc] .= view(MB,1:nc,j)
                end
            #if !flag[j]
                ldiv!(view(MB,:,j),lu(ACCL[ia]),B)
            end
        end
        #ia = findfirst(tM.qwtm.X2x[i,:].*tM.qwtm.T2x[:,t])
    end
end



AD = Matrix(undef,nt,nt);
co = Matrix(undef,nt,nt);
comb0 = CartesianIndex.(zeros(Int64,3),zeros(Int64,3),zeros(Int64,3));
for t1=1:nt
    cf = t1 in findall(flagUpACL)
    ia = findfirst(findlast(flagUpACL[1:t+j-1]).==findall(flagUpACL))
    for t2=1:nt
        #global comb0, comb03
        # if t2==t1
        #     AD[t1,t2] = [cf,0,At[t2],1]
        # end
        # if t2==t1+1
        #     comb = CartesianIndex(1,AD[t1,t2-1][3],At[t2]);
        #     cf2 = comb != comb0
        #     AD[t1,t2] = [cf2,1,At[t2],2]
        #     comb0 = comb
        # end
        # if t2==t1+2
        #     comb3 = CartesianIndex(2,AD[t1,t2-1][3],At[t2]);
        #     cf3 = comb3 != comb03
        #     AD[t1,t2] = [cf3,2,At[t2],3]
        #     comb03 = comb3
        # end
        for j=1:3
            if t1<=t2<=t1+2
                if t2-1 > 0
                    if j-1 >0
                comb = CartesianIndex(j-1,AD[t1,t2-1][3],At[t2]);
                cf2 = comb != comb0[j]
                AD[t1,t2] = [cf2,j-1,At[t2],j-1]
                comb0 = comb[j]
                    end
                end
            end
        end
    end
end


At = zeros(Int64,nt);
ki = findall(flagUpACL)
for (k,v) in enumerate(ki)
    At[v:nt] .= k
end
