n = 10000;
m = 200;
C = zeros(Float64, m)
A = rand(Float64, m,n)
B = rand(Float64, n)


@time for i=1:8000
    mul!(C,A,B)
end

@time for i=1:8000
    C.=A*B
end


@time for i=1:800
    BLAS.gemv('N', A, B)
end
