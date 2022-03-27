ACL = cholesky(mA)
ALU = lu(mA)
ALD = ldlt(mA);

x = similar(b)
xt = [similar(b) for i=1:8]
bt = [copy(b) for i=1:8]
ACLT = [cholesky(mA) for i=1:8]
ALUT = [lu(mA) for i=1:8]
ALDT = [ldlt(mA) for i=1:8]

@time x.=ACL\b
@time x.=ALD\b;

@btime foo($xt,$ACLT,$bt)
@btime fooT($xt,$ACLT,$bt)

@btime foo($xt,$ALUT,$bt)
@btime fooT($xt,$ALUT,$bt)

@btime foo($xt,$ALDT,$bt)
@btime fooT($xt,$ALDT,$bt)

function foo(xt,ACLT,bt)
    for i=1:100
        thr_id = Threads.threadid()
        xt[thr_id] .= ACLT[thr_id]\bt[thr_id];
    end
end

function fooT(xt,ACLT,bt)
    @inbounds Threads.@threads for i=1:100
        thr_id = Threads.threadid()
        xt[thr_id] .= ACLT[thr_id]\bt[thr_id];
    end
end

btt = zeros(length(b),100);
btt.=b
xtt = similar(btt)
@btime $xtt.=$ACL\$btt;

@btime $xtt.=$ALU\$btt;

using SuiteSparse
import SuiteSparse.CHOLMOD: Sparse, Dense, Factor, C_Dense, C_Factor, C_Sparse, SuiteSparse_long
import SuiteSparse.CHOLMOD: common_struct, common_supernodal, common_nmethods, change_factor!, set_print_level, defaults
import SuiteSparse.CHOLMOD: VTypes, @cholmod_name, fact_, cholesky!, CHOLMOD_A



T = Float64
X = [Dense(Array{T}(undef, length(b))) for i in 1:nth]
Eref = [Ref(Ptr{C_Dense{T}}(C_NULL)) for i in 1:nth]
Yref = [Ref(Ptr{C_Dense{T}}(C_NULL)) for i in 1:nth]
B = [Dense(b) for i in 1:nth]

@btime fg($X, $Yref, $Eref, $ACLT, $B)

function fg(X, Yref, Eref, ACLT, B)
  @inbounds Threads.@threads for n=1:100
    thr_id = Threads.threadid()
    sam_solve!(X[thr_id], Yref[thr_id], Eref[thr_id], ACLT[thr_id], B[thr_id])
  end
end

common_struct[3] = fill(0xff, common_size)
common_struct[1][1] = 0x01
common_struct[2][1]

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


AA = make_CL_in_julia(ACL, 4)
@btime zer(xt,AA,bt);
@btime zerT(xt,AA,bt);

function zer(xx,AA,bb)
  for i=1:100
    thr_id = Threads.threadid()
    back_slash_slvr!(xx[thr_id],AA,bb[thr_id])
  end
end

function zerT(xx,AA,bb)
  Threads.@threads for i=1:100
    thr_id = Threads.threadid()
    back_slash_slvr!(xx[thr_id],AA,bb[thr_id])
  end
end
