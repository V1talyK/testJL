using SuiteSparseGraphBLAS
A1 = GBMatrix(A)
b1 = GBMatrix(b)
U = select(TRIU, A1)
L = select(TRIL, A1)
A1*b1
