x0 = zeros(200222)
r0 = b-A*x0
z0 = copy(r0)

@time for i=1:200
    dr0 = dot(r0,r0)
    ak = dr0/dot((A*z0),z0)
    xk = x0+ak*z0
    rk = r0-ak*A*z0
    β = dot(rk,rk)/dr0
    zk = rk+β*z0
    z0[:] = zk
    x0[:] = xk
    r0[:] = rk
    if mod(i,20)==0 println(sum(abs.(b-A*x0))) end
end

x00 = A\b;
mean(abs.(x0-x00))


x0 = zeros(200222)
r0 = b-A*x0
z0 = copy(r0)

@time for i=1:200
    dr0 = dot(r0,r0)
    ak = dr0/dot((A*z0),z0)
    xk = x0+ak*z0
    rk = r0-ak*A*z0
    β = dot(rk,rk)/dr0
    zk = rk+β*z0
    z0[:] = zk
    x0[:] = xk
    r0[:] = rk
    if mod(i,20)==0 println(sum(abs.(b-A*x0))) end
end
