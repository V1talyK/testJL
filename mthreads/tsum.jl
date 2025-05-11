using BenchmarkTools
function sum_single(a)
    s = 0.0
    for i in a
        s += i
    end
    s
end


function sum_multi_good(a)
    chunks = Iterators.partition(a, length(a) ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn sum_single(chunk)
    end
    chunk_sums = fetch.(tasks)
    return sum_single(chunk_sums)
end

ind = rand(10000000);
sum_single(ind)
@btime sum_single($ind);
sum_multi_good(ind)
@btime sum_multi_good($ind)