
function CG_hand!(x,A,b)
    x .= 0
    r0 = A*x.-b
    z0 = copy(r0)
    r = copy(r0);
    z = copy(z0);

    for i=1:1000
        dr = dot(r,r)
        a = dot(r,r)/dot(A*z,z)
        x .= x .+ a*z;
        r .= r .- a*A*z
        bet = dot(r,r)/dr
        z .= r .+ bet*z
    end
end


y = similar(x).*0
CG_hand!(y,A,b)

sum(abs,x.-y)
extrema(x.-y)
