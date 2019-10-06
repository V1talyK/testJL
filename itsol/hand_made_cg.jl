using LinearAlgebra, CuArrays, CuArrays.CUSPARSE
using TimerOutputs
const to = TimerOutput()
reset_timer!(to::TimerOutput)

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

@time cg_cpu(LUi, A,b,1000)
LinearAlgebra.su
function cg_cpu(PCG, A,b,maxIt=1000)
    x = zeros(size(A,1))
    r = b-A*x
    z = copy(r)
    p = copy(z)
    t = copy(z)
    q = copy(z)
    rho = dot(r,z)
    @time for i=1:maxIt
        t[:] = CL.L\r;#LinearAlgebra.BLAS.trsm('R', 'L', 'N', 'N', ones(Float64), CL, b)#CL.L\r
        z[:] = CL.L'\t
        #ldiv!(z,PCG,r)
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
        q[:] = A*p
        temp = dot(p,q)
        α = rho/temp;

        #ak = dr0/dot((A*z0),z0)
        x[:] += α*p
        r[:] -= α*q
        nrmr = norm(r)

        #if mod(i,20)==0 println(sum(abs.(b-A*x))) end
    end
end

cuA = CuArrays.CUSPARSE.CuSparseMatrixCSR(-A[p,p]);
cuB = CuArray(b[p])
CL = cholesky(-A)
CLU = SparseArrays.sparse(CL.L)';
findnz(CL[:U])
cuCLU =  CuArrays.CUSPARSE.CuSparseMatrixCSR(CLU);


infoT = sv_analysis('T', 'T', 'U', icfA,'O')
info = sv_analysis('N', 'T', 'U', icfA,'O')

@time cg_cuda(AU1, cuA,cuB,1000);

function cg_cuda(PCG, cuA,cuB,maxIt=1000)
    infoT = sv_analysis('T', 'T', 'U', PCG,'O')
    info = sv_analysis('N', 'T', 'U', PCG,'O')

    n = length(cuB)
    x = CuArrays.zeros(eltype(cuA),n)
    r = cuB-cuA*x;
    z = copy(r)
    p = copy(z)
    t = copy(z)
    rho = [dot(r,z)]
    q = copy(r)

    cg_cuda_iterator!(x,AU1,r,t,z,p,q,infoT,info,maxIt,rho,n)
    return x
end

function cg_cuda_iterator!(x,AU1,r,t,z,p,q,infoT,info,maxIt,rho,n)
    @time @inbounds for i=1:maxIt
        #@timeit to "t" t[:] .= CUSPARSE.sv_solve('T','U',one(Float64),AU1,r,infoT,'O')
        @timeit to "t1"  CUSPARSE.sv_solve!('T','U',one(Float64),AU1,r,t,infoT,'O')
        #r1=copy(r)
        #@timeit to "t2"  CUSPARSE.sv2!('T','U',one(Float64),AU1,r1,'O')
        # r_cpu = Array(r)
        # @timeit to "t_cpu" AU\r_cpu
        @timeit to "z" CUSPARSE.sv_solve!('N','U',one(Float64),AU1,t,z,info,'O')
        #t = cuCLL\r;
        #z = cuCLL'\t
        #global rho
        #rho = dot(r,r)
        @timeit to "t2" rhop = copy(rho)
        @timeit to "t3" rho[1] = CuArrays.CUBLAS.dot(n,r,1,z,1);
        @timeit to "if"         if i!=1
            β = rho[1]/rhop[1]
            #z[:] = z+β*p
            CuArrays.CUBLAS.axpy!(n,β,p,1,z,1) #z[:] = z+β*p
        end
        @timeit to "t4" p[:] = copy(z)
        #CuArrays.CUBLAS.cublasDcopy_v2(CuArrays.CUBLAS.handle(),n,r,1,z,1)
        #y = α ∗ op ( A ) ∗ x + β ∗ y
        #q[:] = cuA*p
        @timeit to "mv"         CuArrays.CUSPARSE.mv!('N',one(Float64),cuA,p,zero(Float64),q,'O')
        #temp = dot(p,q)
        @timeit to "t5" temp = CuArrays.CUBLAS.dot(n,p,1,q,1);
        @timeit to "t6" α = rho[1]/temp;

        #ak = dr0/dot((A*z0),z0)
        @timeit to "t7" CuArrays.CUBLAS.axpy!(n,α,p,1,x,1) #x[:] += α*p
        @timeit to "t8" CuArrays.CUBLAS.axpy!(n,-α,q,1,r,1) #r[:] += -α*q
        nrmr = norm(r)
        #println(nrmr)
        #if mod(i,100)==0 println(sum(abs.(cuB-cuA*x))) end
    end
end
CuArrays.CUBLAS.axpy!()
