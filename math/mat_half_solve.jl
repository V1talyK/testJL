A = rand(3,3)
B = rand(3,2)
C = rand(2,3)
D = rand(2,2)
B1 = ones(3)
B2 = 1 .+ones(2)

XY0 = [A B;C D]\[B1; B2]


(D\B')'*B2
(A .- (D'\B')'*C)\(B1 .- (D'\B')'*B2 )
