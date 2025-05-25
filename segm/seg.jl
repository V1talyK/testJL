using Clustering, Plots

nt  = 100


using TSML

function top_down_segmentation(ts::Vector{Float64}, max_error::Float64, min_length::Int)
    segments = []
    recursive_split!(ts, 1, length(ts), max_error, min_length, segments)
    return segments
end

function recursive_split!(ts, start_idx, end_idx, max_error, min_length, segments)
    if end_idx - start_idx + 1 < 2*min_length
        push!(segments, (start_idx, end_idx))
        return
    end
    
    best_split = find_best_split(ts, start_idx, end_idx)
    error_left = calculate_approximation_error(ts, (start_idx, best_split))
    error_right = calculate_approximation_error(ts, (best_split+1, end_idx))
    
    if error_left + error_right > max_error
        recursive_split!(ts, start_idx, best_split, max_error, min_length, segments)
        recursive_split!(ts, best_split+1, end_idx, max_error, min_length, segments)
    else
        push!(segments, (start_idx, end_idx))
    end
end

function find_best_split(ts, start_idx, end_idx)
    min_error = Inf
    best_split = start_idx
    
    for split_point in start_idx+1:end_idx-1
        error_left = calculate_approximation_error(ts, (start_idx, split_point))
        error_right = calculate_approximation_error(ts, (split_point+1, end_idx))
        total_error = error_left + error_right
        
        if total_error < min_error
            min_error = total_error
            best_split = split_point
        end
    end
    
    return best_split
end

function calculate_approximation_error(ts, segment)
    start_idx, end_idx = segment
    segment_data = ts[start_idx:end_idx]
    # Simple linear approximation
    x = 1:length(segment_data)
    coeff = linreg(x, segment_data)
    approx = coeff[1] .+ coeff[2] .* x
    return sum((segment_data .- approx).^2)
end

function linreg(x, y)
    n = length(x)
    A = [ones(n) x]
    println("--")
    display(A'*A)
    println(y)
    return (A' * A) \ (A' * y)
end

nt = 256
vt = 1:nt
ts = vcat(0.1*vt[1:128],0.1*vt[129].+0.2*vt[1:128])
plot(ts)

 se = top_down_segmentation(ts, 1e-4, 2)