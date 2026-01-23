using SparseArrays, LinearAlgebra
using Plots

A = sprand(100,100,0.1)
A = A'*A
Ac = sum(A, dims = 1)[:]
for i = 1:100
    A[i, i] = -Ac[i]
end
iA = inv(Matrix(A))
U, E, V = svd(iA)

iAr = U * Diagonal(E) * V'
nn = 50
iAr2 = U[:,1:nn]*Diagonal(E[1:nn])*V'[1:nn,:]

sum(abs2, iAr.-iA)
sum(abs2, iAr2.-iA)

bb = rand(100)
sum(abs2, iA*bb .- iAr2*bb)
scatter(iA*bb, iAr2*bb)

0.18*0.024*2

1.350*10/(0.18*0.024*2)

1350/(8*2.8*2)