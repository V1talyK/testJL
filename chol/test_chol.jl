using LinearAlgebra, BenchmarkTools, Profile
n = 100
A = rand(n,n); A = A.+A'; A = diagm(sum(A,dims=2)[:].+1).-A
mA = sprand(n,n,0.3); mA = mA.+mA'; mA = diagm(sum(mA,dims=2)[:].+1).-mA
b = rand(n)


function hand_chol(A)
    L = zeros(size(A))
    n = size(A,1)
    for j = 1:n
        L[j,j] = sqrt(A[j,j] - sum(L[j,1:j-1].^2))
        for i = j+1:n
            L[i,j] = 1/L[j,j]*(A[i,j]-sum(L[i,1:j-1].*L[j,1:j-1]))
        end
    end
    return L
end

@time L = hand_chol(A)
@time cholesky(A)
sum(abs, L*L'.-A)


@btime hand_cholS($mA);
@btime hand_chol($mA);
@time cholesky(mA);
L = hand_cholS(mA);
sum(abs, L*L'.-mA)

A = copy(mA)
