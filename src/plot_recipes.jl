@recipe function f(d::UvModelDensity, interval::Tuple, npoints::Int = 100)
    xs = range(interval[1], stop = interval[2], length = npoints)
    ys = map(x -> evaluate(d, x), xs) 
    xs, ys
end
@recipe function f(d::UvModelDensity, interval::Interval, npoints::Int = 100)
    xs = range(interval.left, stop = interval.right, length = npoints)
    ys = map(x -> evaluate(d, x), xs) 
    xs, ys
end
@recipe function f(d::UvModelDensity, h::Histogram{<:Any,1}, npoints::Int = 5*length(h.edges[1]))
    @assert isa(h.edges[1], AbstractRange) "Plotting for uneven spaced histgrams is not yet implemented. "
    xs = range(h.edges[1][1], stop = h.edges[1][end], length = npoints)
    ys = map(x -> evaluate(d, x), xs) 
    xs, ys .* step(h.edges[1])
end