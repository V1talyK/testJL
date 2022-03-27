import SuiteSparse.CHOLMOD: Sparse, Dense, Factor, C_Dense, C_Factor, C_Sparse, SuiteSparse_long
import SuiteSparse.CHOLMOD: common_struct, common_supernodal, common_nmethods, change_factor!, set_print_level, defaults
import SuiteSparse.CHOLMOD: VTypes, @cholmod_name, fact_, cholesky!, CHOLMOD_A
ACL = cholesky(mA)
@time x .= ACL\b;

nth = Threads.nthreads()
T = Float64

X = [Dense(zeros(T, length(b))) for i in 1:nth]
Eref = [Ref(Ptr{C_Dense{T}}(C_NULL)) for i in 1:nth]
Yref = [Ref(Ptr{C_Dense{T}}(C_NULL)) for i in 1:nth]

bd = [Dense(b) for i in 1:nth];
#Threads.@threads
ACLT = [cholesky(mA) for i in 1:nth]
@time @sync for i=1:100
  Threads.@spawn bar(X,Yref,Eref,ACLT,bd)
end

@time  bar(100,X,Yref,Eref,ACLT,bd)
@time  barT(100,X,Yref,Eref,ACLT,bd)
@btime bar(100,$X,$Yref,$Eref,$ACLT,$bd)
@btime barT(100,$X,$Yref,$Eref,$ACLT,$bd)


bb = [copy(b) for i in 1:nth]
xx = [copy(b) for i in 1:nth]

@btime foo($200,$xx,$ACLT,$bb)
@btime fooT($200,$xx,$ACLT,$bb)

function foo(n,x,ACLT,b)
  for i=1:n
      thr_id = Threads.threadid()
      x[thr_id] .= ACLT[thr_id]\b[thr_id]
  end
end

function fooT(n,x,ACLT,b)
  Threads.@threads for i=1:n
    #Threads.@spawn begin
      thr_id = Threads.threadid()
      x[thr_id] .= ACLT[thr_id]\b[thr_id]
    #end
  end
end

function bar(n, X,Yref,Eref,ACLT,bd)
  for i=1:n
      thr_id = Threads.threadid()
      sam_solve!(X[thr_id], Yref[thr_id], Eref[thr_id], ACLT[thr_id], bd[thr_id]);
  end
  return nothing
end


function barT(n, X,Yref,Eref,ACLT,bd)
  @sync for i=1:n
      Threads.@spawn bar2(X,Yref,Eref,ACLT,bd)
  end
  return nothing
end


function bar2(X,Yref,Eref,ACLT,bd)
  thr_id = Threads.threadid()
  sam_solve!(X[thr_id], Yref[thr_id], Eref[thr_id], ACLT[thr_id], bd[thr_id]);
end


function sam_solve!(X::Dense{Tv}, Yref::Ref{Ptr{C_Dense{Tv}}}, Eref::Ref{Ptr{C_Dense{Tv}}}, F::Factor{Tv}, B::Dense{Tv}) where Tv<:VTypes
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

function copyto!(dest::Dense{T},src::Vector{T}) where {T<:AbstractFloat}
  #T = eltype(dest)
  GC.@preserve dest begin
    s = unsafe_load(pointer(dest));
    pt = Ptr{T}(s.x)
    @inbounds for (i, c) in enumerate(eachindex(src))
        #@time x[c];
        unsafe_store!(pt, src[c], i);
    end
  end
end
