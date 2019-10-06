Ao  = A.rowval;
Ac  = A.colptr;
Av  = A.nzval;

function color_jpl_kernel(n::Int64, c::Int64, Ao::Array{<:Signed,1},
                                 Ac::Array{<:Signed,1}, Av::Array{<:AbstractFloat,1},
                                 randoms::Array{<:Signed,1}, colors::Array{<:Signed,1})

  for i=1:n
    f=true; # true iff you have max random
    # ignore nodes colored earlier
    if ((colors[i] != 0)) continue; end;

    ir = randoms[i];

    # look at neighbors to check their random number
    for k = Ao[i]:Ao[i+1]-1
      # ignore nodes colored earlier (and yourself)
      j = Ac[k];
      jc = colors[j];
      if (((jc != -1) && (jc != c)) || (i == j))
        continue;
      end

      jr = randoms[j];
      if ir <= jr         f=false;          end

    end

    # assign color if you have the maximum random number
    if f  colors[i] = c;    end;

  end
end
#define CUDA_MAX_BLOCKS <Maximum blocks to launch, depends on GPU>

function color_jpl(n::Int64,Ao::Array{<:Signed,1},
                  Ac::Array{<:Signed,1}, Av::Array{<:AbstractFloat,1})

    randoms = rand(1:3,n); # allocate and init random array
    colors = fill(0, n);  # init colors to -1

    for c=1:n
        color_jpl_kernel(n, c, Ao, Ac, Av,randoms, colors)
        left = count(i->i==0,colors);
        if left == 0 break; end;
    end
    return colors
end

n = size(A,1)
clr = color_jpl(n,Ao,Av,Av)

cuA =  CuArrays.CUSPARSE.CuSparseMatrixCSR(A)
Av = Array(cuA.nzVal)
Ao = Array(cuA.rowPtr)
Ac = Array(cuA.colVal)

p = sortperm(clr)
