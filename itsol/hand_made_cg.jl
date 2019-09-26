using LinearAlgebra, CuArrays, CuArrays.CUSPARSE
x0 = zeros(size(A,1))
@time r0 = b-A*x0
Y = copy(b)
@time BLAS.sbmv!(uplo, k, 1, A, x0, 1, Y)
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

@time x00 = A\b;
mean(abs.(x0-x00))

A = -A
CL = cholesky(A)
x = zeros(size(A,1))
r = b-A*x
z = copy(r)
p = copy(z)
rho = dot(r,z)
@time for i=1:200
    t = CL.L\r
    z = CL.L'\t
    global rho
    #rho = dot(r,r)
    rhop = copy(rho)
    rho = dot(r,z)
    if i==1
        p[:] = copy(z)
    else
        β = rho/rhop
        z[:] = z+β*p
        p[:] = z
    end


    #y = α ∗ op ( A ) ∗ x + β ∗ y
    q = A*p
    temp = dot(p,q)
    α = rho/temp;

    #ak = dr0/dot((A*z0),z0)
    x[:] += α*p
    r[:] -= α*q
    nrmr = norm(r)

    if mod(i,20)==0 println(sum(abs.(b-A*x))) end
end

cuA = CuArrays.CUSPARSE.CuSparseMatrixCSR(-A);
cuB = CuArray(b)
CL = cholesky(-A)
CLU = SparseArrays.sparse(CL.L)';
findnz(CL[:U])
cuCLU =  CuArrays.CUSPARSE.CuSparseMatrixCSR(CLU);

infoT = sv_analysis('T', 'T', 'U', AU1,'O')
infoT = sv_analysis('T', 'T', 'U', icfA,'O')
info = sv_analysis('N', 'T', 'U', AU1,'O')
info = sv_analysis('N', 'T', 'U', icfA,'O')

x = CuArrays.zeros(eltype(cuA),size(A,1))
r = cuB-cuA*x;
r=r
z = copy(r)
p = copy(z)
rho = [dot(r,z)]
q = copy(r)

for i=1:5000
    t = CUSPARSE.sv_solve('T','U',one(Float64),AU1,r,infoT,'O')
    z[:] = CUSPARSE.sv_solve('N','U',one(Float64),AU1,t,info,'O')
    #t = cuCLL\r;
    #z = cuCLL'\t
    #global rho
    #rho = dot(r,r)
    rhop = copy(rho)
    rho[1] = dot(r,z)
    if i==1
        p[:] = copy(z)
    else
        β = rho[1]/rhop[1]
        #z[:] = z+β*p
        z[:] = z+β*p
        p[:] = copy(z)
    end

    #y = α ∗ op ( A ) ∗ x + β ∗ y
    q[:] = cuA*p
    temp = dot(p,q)
    α = rho[1]/temp;

    #ak = dr0/dot((A*z0),z0)
    x[:] += α*p
    r[:] += -α*q
    nrmr = norm(r)

    if mod(i,20)==0 println(sum(abs.(cuB-cuA*x))) end
end
