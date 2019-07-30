function mean(h::Histogram{T1, 1}) where {T1 <: Real}
    mean(midpoints(h.edges[1]), FrequencyWeights(h.weights))
end

function std(h::Histogram{T1, 1}, μ = mean(h)) where {T1 <: Real}
    std(midpoints(h.edges[1]), FrequencyWeights(h.weights), mean = μ, corrected = true)
end

function mean_and_std(h::Histogram{T1, 1}) where {T1 <: Real}
    μ = mean(h)
    σ = std(h, μ = μ)
    return μ, σ 
end