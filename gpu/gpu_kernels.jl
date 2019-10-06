function dot1(a,b,c, N, threadsPerBlock, blocksPerGrid)

    # Set up shared memory cache for this current block.
    cache = @cuDynamicSharedMem(Float64, threadsPerBlock)

    # Initialise some variables.
    tid = (threadIdx().x - 1) + (blockIdx().x - 1) * blockDim().x
    totalThreads = blockDim().x * gridDim().x
    cacheIndex = threadIdx().x - 1
    temp = 0

    # Iterate over vector to do dot product in parallel way
    while tid < N
        temp += a[tid + 1] * b[tid + 1]
        tid += totalThreads
    end

    # set cache values
    cache[cacheIndex + 1] = temp

    # synchronise threads
    sync_threads()

    # In the step below, we add up all of the values stored in the cache
    i::Int = blockDim().x ÷ 2
    while i!=0
        if cacheIndex < i
            cache[cacheIndex + 1] += cache[cacheIndex + i + 1]
        end
        sync_threads()
        i = i ÷ 2
    end

    # cache[1] now contains the sum of vector dot product calculations done in
    # this block, so we write it to c
    if cacheIndex == 0
        c[blockIdx().x] = cache[1]
    end

    return nothing
end

function sum1(iter_result, iter_num)
    if threadIdx.x == 0
        atomicAdd(iter_result + iter_num, block_ressult);
    end
    if (cacheIndex==0) {
    double result=cuda_atomicAdd(dot_res,cache[0]);
}
    return nothing
end

function main()

    # Initialise variables
    N::Int64 = 33 * 1024
    threadsPerBlock::Int64 = 256
    blocksPerGrid::Int64 = min(32, (N + threadsPerBlock - 1) / threadsPerBlock)

    # Create a,b and c
    a = CuArrays.CuArray(fill(0, N))
    b = CuArrays.CuArray(fill(0, N))
    c = CuArrays.CuArray(fill(0, blocksPerGrid))

    # Fill a and b
    for i in 1:N
        a[i] = i
        b[i] = 2*i
    end

    # Execute the kernel. Note the shmem argument - this is necessary to allocate
    # space for the cache we allocate on the gpu with @cuDynamicSharedMem
    @cuda blocks = blocksPerGrid threads = threadsPerBlock shmem =
    (threadsPerBlock * sizeof(Int64)) dot1(a,b,c, N, threadsPerBlock, blocksPerGrid)

    # Copy c back from the gpu (device) to the host
    c = Array(c)

    local result = 0

    # Sum the values in c
    for i in 1:blocksPerGrid
        result += c[i]
    end

    # Check whether output is correct
    println("Does GPU value ", result, " = ", 2 * sum_squares(N - 1))
end
function main1()

    # Initialise variables
    N::Int64 = 1e6
    threadsPerBlock::Int64 = 256
    blocksPerGrid::Int64 = min(32, (N + threadsPerBlock - 1) / threadsPerBlock)

    # Create a,b and c
    a = CuArray(rand(Float64,N))
    b = CuArray(rand(Float64,N))
    c = CuArrays.CuArray(fill(0., blocksPerGrid))

    # Execute the kernel. Note the shmem argument - this is necessary to allocate
    # space for the cache we allocate on the gpu with @cuDynamicSharedMem
    @time begin
        @cuda blocks = blocksPerGrid threads = threadsPerBlock shmem =
        (threadsPerBlock * sizeof(Int64)) dot1(a,b,c, N, threadsPerBlock, blocksPerGrid)
    end
        # Copy c back from the gpu (device) to the host
        c = Array(c)

        local result = 0

        # Sum the values in c
        for i in 1:blocksPerGrid
            result += c[i]
        end

    # Check whether output is correct
    println("Does GPU value ", result, " = ", dot(a,b))
end

function sum_squares(x)
    return (x * (x + 1) * (2 * x + 1) / 6)
end

main()
main1()

N = Int64(1e6)
threadsPerBlock = 256
blocksPerGrid  = Int64(min(32, (N + threadsPerBlock - 1) / threadsPerBlock))

r = CuArray(rand(Float64,N))
z = CuArray(rand(Float64,N))
rsl = CuArray(fill(0,blocksPerGrid))
@cuda blocks = blocksPerGrid threads = threadsPerBlock shmem =(threadsPerBlock * sizeof(Int64)) dot1(r,z,rsl, N, threadsPerBlock, blocksPerGrid)

@cuda threads=12 sum1(r, rt)
