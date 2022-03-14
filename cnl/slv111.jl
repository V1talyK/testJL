using SparseArrays
using LinearAlgebra
using SuiteSparse
using BenchmarkTools

import SuiteSparse.CHOLMOD: Sparse, Dense, Factor, C_Dense, C_Factor, C_Sparse, SuiteSparse_long
import SuiteSparse.CHOLMOD: common_struct, common_supernodal, common_nmethods, change_factor!, set_print_level, defaults
import SuiteSparse.CHOLMOD: VTypes, @cholmod_name, fact_, cholesky!, CHOLMOD_A
import SparseArrays: getcolptr, getrowval

function supercholesky(A::Sparse; shift::Real=0.0, check::Bool = true,
  perm::Union{Nothing,AbstractVector{SuiteSparse_long}}=nothing)

  cm = defaults(common_struct[Threads.threadid()])
  set_print_level(cm, 0)

  # Force a supernodal solution (eliminates alloc on solve)
  unsafe_store!(common_supernodal[Threads.threadid()], 2)

  # Compute the symbolic factorization
  F = fact_(A, cm; perm = perm)

  # Compute the numerical factorization
  cholesky!(F, A; shift = shift, check = check)

  return F
end




T = Float64
N = length(b)

X = Dense(Array{T}(undef, N))
B = Dense(b)
Yref = Ref(Ptr{C_Dense{T}}(C_NULL))
Eref = Ref(Ptr{C_Dense{T}}(C_NULL))

As = Sparse(A)

F_simple = cholesky(A)
F_super =  supercholesky(As);

_solve!(X, Yref, Eref, F_simple, B);
@show F_simple
A\B ≈ X ? println("Simplical pass") : error("Simplical fail")

_solve!(X, Yref, Eref, F_super, B);
@show F_super
A\B ≈ X ? println("Supernodal pass") : error("Supernodal fail")

b1 = @benchmark _solve!($X, $Yref, $Eref, $F_simple, $B)

b2 = @benchmark _solve!($X, $Yref, $Eref, $F_super, $B)

b3 = @benchmark _mysolve!($X, $Yref, $Eref, $F_simple, $B)

b4 = @benchmark _mysolve!($X, $Yref, $Eref, $F_super, $B)



function _mysolve!(X::Dense{Tv}, Yref::Ref{Ptr{C_Dense{Tv}}}, Eref::Ref{Ptr{C_Dense{Tv}}}, F::Factor{Tv}, B::Dense{Tv}) where Tv<:VTypes
  Y = unsafe_load(Yref[])
  nrow = Int(Y.nzmax / Y.ncol)
  d = nrow
  unsafe_store!(Yref[], C_Dense{Float64}(nrow, Y.ncol, Y.nzmax, d, Y.x, Y.z, Y.xtype, Y.dtype))
  # Pointer to pre-allocated dense matrix
  Xref = Ptr{C_Dense{Tv}}(pointer(X))
  sys = CHOLMOD_A # Solve type
  Bset = C_NULL   # Unused parameter
  Xset = C_NULL   # Unused parameter

  if size(F,1) != size(B,1)
    throw(DimensionMismatch("LHS and RHS should have the same number of rows. " *
    "LHS has $(size(F,1)) rows, but RHS has $(size(B,1)) rows."))
  end

  if !issuccess(F)
    s = unsafe_load(pointer(F))
    if s.is_ll == 1
      throw(LinearAlgebra.PosDefException(s.minor))
    else
      throw(LinearAlgebra.ZeroPivotException(s.minor))
    end
  end

  res = ccall((@cholmod_name("solve2"), :libcholmod), Cint,
  (Cint, Ptr{C_Factor{Tv}}, Ptr{C_Dense{Tv}}, Ptr{C_Sparse{Tv}}, Ref{Ptr{C_Dense{Tv}}},  Ref{Ptr{C_Sparse{Tv}}},  Ref{Ptr{C_Dense{Tv}}},  Ref{Ptr{C_Dense{Tv}}}, Ptr{UInt8}),
  sys,   F, B,Bset,Xref,Xset,Yref,Eref,SuiteSparse.CHOLMOD.common_struct[Threads.threadid()])

  if (res != 1)
    throw(ErrorException("CHOLMOD solve failure"))
  end

  return nothing
end


@time ACL = cholesky(mA);
@time cholesky!(ACL,mA);
x = similar(b)
@time x .= ACL\b;

@time sA = Sparse(mA);
@time sACL = supercholesky(sA)
@time cholesky!(sACL,sA);

@time F = fact_(sA, cm; perm = perm);
@time cholesky!(ACL, sA)


x = zeros(10)
xd = Dense(x)
xd[1]

pp = pointer(xd)
unsafe_store!(pp, 2.,1)
