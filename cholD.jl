CA = cholesky(-AS)
SuiteSparse.CHOLMOD.lowrankupdate(CA, v1)

using SuiteSparse

ACL0 = preInterpPPw(nc,grdD,GM[:,1],GG,AG,actW[:,t],WIv,bm,BGQ,cB,da[t])
WIv66 = WIv.*wst;
t = 68
@time ACL68 = preInterpPPw(nc,grdD,GM[:,1],GG,AG,actW[:,t],WIv,bm,BGQ,cB,da[t])
WIv68 = WIv.*wst;

WIv = WI[tt][actW[:,66]];
AW66 = makeAWellChol(actW[:,66],grdD["Won"],WIv,nc)

WIv = WI[tt][actW[:,68]];
AW68 = makeAWellChol(actW[:,t],grdD["Won"],WIv,nc)


sum(AW68 - AW66)
AW66m = zeros(2510);
@time AW66m = sqrt.(AW66);
@time ACL066 = SuiteSparse.CHOLMOD.lowrankupdate(ACL0, AW66m)
sum(ACL66\ones(2510)) - sum(ACL066\ones(2510))

ACL066 = copy(ACL0)

@time udateCL!(ACL066,AW66m)

function udateCL!(CL,v0)
    v = zeros(eltype(v0),length(v0));
    ia = findall(v0.!=0)
    for i in ia
        v[i] = v0[i]
        SuiteSparse.CHOLMOD.lowrankupdate!(CL, v)
        v[i]=0;
    end
end
