# This file is a part of RadiationSpectra.jl, licensed under the MIT License (MIT).

using Test
using Random
using RadiationSpectra
using RadiationSpectra: match_calibration_lines, determine_calibration_constant_through_peak_ratios

# Exhaustive search over monotone one-to-one matchings: most matches, then smallest sum of squared relative residuals.
function matching_oracle(positions, lines, c; rtol, min_match_span)
    n, m = length(positions), length(lines)
    min_span = min_match_span * (lines[end] - lines[1])
    best = (0, -Inf)
    for peak_set in 1:(2^n - 1), line_set in 1:(2^m - 1)
        count_ones(peak_set) == count_ones(line_set) || continue
        peaks = [i for i in 1:n if (peak_set >> (i - 1)) & 1 == 1]
        matched = [j for j in 1:m if (line_set >> (j - 1)) & 1 == 1]
        lines[matched[end]] - lines[matched[1]] >= min_span || continue
        residuals = [(c * positions[i] - lines[j]) / lines[j] for (i, j) in zip(peaks, matched)]
        all(abs.(residuals) .<= rtol) || continue
        best = max(best, (length(peaks), -sum(abs2, residuals)))
    end
    return (best[1], -best[2])
end

@testset "Auto calibration" begin
    @testset "Line matching" begin
        rng = Xoshiro(42)
        for _ in 1:100
            lines = sort(100 .+ 900 .* rand(rng, rand(rng, 3:5)))
            c = exp(5 * rand(rng) - 4)
            rtol = rand(rng, (0.005, 0.01, 0.03))
            present = filter(_ -> rand(rng, Bool), lines)
            noise = 1 .+ 1.6 * rtol .* (2 .* rand(rng, length(present)) .- 1)
            spurious = 100 .+ 900 .* rand(rng, rand(rng, 0:3))
            positions = sort(unique(vcat(present .* noise, spurious) ./ c))
            min_match_span = rand(rng, (0.0, 0.3, 0.5, 0.8))
            matches = match_calibration_lines(positions, lines, c; rtol, min_match_span)
            peaks, matched = matches.peak_indices, matches.line_indices
            residuals = [(c * positions[i] - lines[j]) / lines[j] for (i, j) in zip(peaks, matched)]
            @test length(peaks) == length(matched) && issorted(peaks, lt=<=) && issorted(matched, lt=<=)
            @test all(abs.(residuals) .<= rtol + 8eps(Float64))
            best_n, best_cost = matching_oracle(positions, lines, c; rtol, min_match_span)
            @test length(peaks) == best_n
            if isempty(peaks)
                @test matches.score == Inf
            else
                @test lines[matched[end]] - lines[matched[1]] >= min_match_span * (lines[end] - lines[1])
                @test sum(abs2, residuals) ≈ best_cost
                @test matches.score ≈ best_cost / best_n
            end
        end

        # An assignment that meets the span requirement wins over a narrower one that fits better.
        @test determine_calibration_constant_through_peak_ratios([100.49, 110.539, 149.4], [100.0, 110.0, 149.0, 150.0, 200.0]; rtol=0.005, min_n_matches=3, min_match_span=0.5) ≈ 1.0 rtol=0.005
        # Equally spaced lines leave the gain ambiguous.
        @test_throws "Ambiguous calibration" determine_calibration_constant_through_peak_ratios([10, 20, 30, 40, 100, 200, 300, 400.0], [1.0, 2.0, 3.0, 4.0]; min_n_matches=4)
        @test_throws ArgumentError determine_calibration_constant_through_peak_ratios([100.0, 200.0, 300.0], [1.0, 2.0, 3.0]; rtol=2)
    end
end
