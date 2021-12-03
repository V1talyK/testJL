using SuiteSparseGraphBLAS
A1 = GBMatrix(A)
b1 = GBMatrix(b)
U = select(SelectOps.TRIU, A1)
L = select(SelectOps.TRIL, A1)
A1\b1

x2 = IterativeSolvers.cg(A1, b1);

function cohen(A)
         U = select(triu, A)
         L = select(tril, A)
         return reduce(+, mul(L, U, (+, pair); mask=A)) ÷ 2
       end

function sandia(A)
  L = select(tril, A)
  return reduce(+, mul(L, L, (+, pair); mask=L))
end

A = GBMatrix([1,1,2,2,3,4,4,5,6,7,7,7], [2,4,5,7,6,1,3,6,3,3,4,5], [1:12...])
M = eadd(A, A', +)
 cohen(A)
12
