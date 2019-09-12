function csrlsvlu!(A::CuSparseMatrixCSR{$elty},
                   b::CuVector{$elty},
                   x::CuVector{$elty},
                   tol::$relty,
                   reorder::Cint,
                   inda::Char)
    cuinda = cusparseindex(inda)
    n = size(A,1)
    if size(A,2) != n
        throw(DimensionMismatch("LU factorization is only possible for square matrices!"))
    end
    if size(A,2) != length(b)
        throw(DimensionMismatch("second dimension of A, $(size(A,2)), must match the length of b, $(length(b))"))
    end
    if length(x) != length(b)
        throw(DimensionMismatch("length of x, $(length(x)), must match the length of b, $(length(b))"))
    end

    #Mat     = similar(A)
    #transpose!(Mat, A)
    cudesca = cusparseMatDescr_t(CUSPARSE_MATRIX_TYPE_GENERAL, CUSPARSE_FILL_MODE_LOWER, CUSPARSE_DIAG_TYPE_NON_UNIT, cuinda)
    rcudesca = Ref{cusparseMatDescr_t}(cudesca)
    singularity = Ref{Cint}(1)
    @check ccall(($(string(fname)),libcusolver), cusolverStatus_t,
                      (cusolverSpHandle_t, Cint, Cint,
                       Ptr{cusparseMatDescr_t}, CuPtr{$elty}, CuPtr{Cint},
                       CuPtr{Cint}, CuPtr{$elty}, $relty, Cint, CuPtr{$elty},
                       Ptr{Cint}),
                      sparse_handle(), n, A.nnz, rcudesca,
                      A.nzVal, A.rowPtr, A.colVal, b, tol, reorder,
                      x, singularity)

    if singularity[] != -1
        throw(SingularException(singularity[]))
    end

    x
end
