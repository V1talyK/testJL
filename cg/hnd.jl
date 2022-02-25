using LinearAlgebra

function CG_hand!(x,A,b)
    x .= 0
    r0 = A*x.-b
    r = copy(r0);
    z = copy(r0);
    tmp = similar(x)
    for i=1:200
        dr = dot(r,r)
        mul!(tmp,A,z)
        a = dr/dot(tmp,z)
        x .= x .+ a.*z;
        r .= r .- a.*tmp
        bet = dot(r,r)/dr
        z .= r .+ bet.*z
    end
end


y = similar(x).*0
@btime CG_hand!(y,A,b)
@profiler CG_hand!(y,A,b)

sum(abs,x.+y)
extrema(x.+y)
ACL = cholesky(mA)

function CG_hand!(x,A,b)
    x .= 0
    r = b - mA*x
    mul!(z, mA,r)
    p = copy(z)
    tmp = similar(x)
    for i=1:200
        a = dot(r,z)
        x .= x + a*p
        r = r + a*A*p
        z = ACL\r

        dr = dot(r,r)
        mul!(tmp,A,z)
        a = dr/dot(tmp,z)
        x .= x .+ a.*z;
        r .= r .- a.*tmp
        bet = dot(r,r)/dr
        z .= r .+ bet.*z
    end
end
