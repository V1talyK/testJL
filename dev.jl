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
