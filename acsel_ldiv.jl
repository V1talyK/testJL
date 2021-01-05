@time x0 = A\b;
ACL = cholesky(mA)
@time ACL\b

L = sparse(ACL.L)
@time Lt = LowerTriangular(L)
U = copy(L');
@time Ut = UpperTriangular(U)

bp = b[ACL.p];

@time y = Lt\bp
@time x .= Ut\y

@time begin
    ldiv!(y,Lt,bp)
    ldiv!(x,Ut,y)
end

sum(abs,x.+x0[ACL.p])
@time L\b;
