# This file is a part of RadiationSpectra.jl, licensed under the MIT License (MIT).

using Test
using Random, StatsBase
using RadiationSpectra
using RadiationSpectra: match_calibration_lines, select_calibration_solution, determine_calibration_constant_through_peak_ratios, calibrate_spectrum, peakfinder

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

# Gaussian peaks of the given amplitudes on an exponential background, rounded to counts.
function peak_histogram(edges, peaks; σ, background=20.0)
    weights = [background * exp(-x / last(edges)) + sum(a * exp(-(x - μ)^2 / (2σ^2)) for (μ, a) in peaks; init=0.0) for x in StatsBase.midpoints(edges)]
    return Histogram(edges, round.(weights), :left)
end

peakfinder_result(h; kwargs...) = (pf = peakfinder(h; kwargs...); (pf[1].weights, pf[2]))

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

    @testset "Peak ratios on histograms" begin
        h = peak_histogram(0:1.0:300, [(100.0, 500.0), (200.0, 500.0)]; σ=3.0, background=2.0)
        c, h_deconv, positions, threshold = determine_calibration_constant_through_peak_ratios(h, [1.0, 2.0]; σ=3.0)
        @test c ≈ 0.01 rtol=0.01
        @test (h_deconv.weights, positions) == peakfinder_result(h; threshold, σ=3.0)
        @test_throws "No peaks identified" determine_calibration_constant_through_peak_ratios(h, [1.0, 2.0]; σ=3.0, min_n_peaks=10)
        @test_throws ArgumentError determine_calibration_constant_through_peak_ratios(h, [1.0, 2.0]; σ=3.0, min_n_peaks=10, rtol=2)
        @test_throws ArgumentError determine_calibration_constant_through_peak_ratios(h, [1.0, 2.0]; σ=3.0, threshold=Inf)

        # Two-match searches return all raw positions, including those the shape filter or max_n_peaks drop,
        # and all peaks of the selected threshold compete for the gain.
        h_pe = peak_histogram(0:0.2:30, [(4.0, 100.0), (8.0, 300.0), (12.0, 100.0), (16.0, 300.0), (20.0, 100.0)]; σ=0.6, background=5.0)
        c, h_deconv, positions, threshold = determine_calibration_constant_through_peak_ratios(h_pe, collect(1.0:5); min_n_peaks=2, max_n_peaks=3, σ=3.0, threshold=10, α=0.05, rtol=0.1)
        @test c ≈ 0.25 rtol=0.01
        @test length(positions) == 5
        @test (h_deconv.weights, positions) == peakfinder_result(h_pe; threshold, σ=3.0)

        # Equal support at different thresholds is not ambiguous, the higher threshold wins.
        solution(c, threshold) = (; c, n_matches=3, score=1e-6, threshold)
        @test select_calibration_solution([solution(1.0, 7.5), solution(2.0, 10.0)]; rtol=0.005, min_n_matches=3).c == 2.0
        @test_throws "Ambiguous calibration" select_calibration_solution([solution(1.0, 10.0), solution(2.0, 10.0)]; rtol=0.005, min_n_matches=3)
    end

    @testset "Spectrum calibration" begin
        lines = [609.312, 911.204, 1120.287, 1460.830, 1764.494, 2614.533]
        # A line missing from the spectrum must neither be matched to a nearby peak nor be fitted.
        for (present, decoys) in ((lines, Float64[]), (lines[[1, 2, 4, 5, 6]], [1120.287 + 15]))
            h = peak_histogram(0:50.0:280_000, [(l / 0.01, 500.0) for l in vcat(present, decoys)]; σ=170.0)
            _, h_deconv, positions, threshold, c, c_precal = calibrate_spectrum(h, lines; σ=3.0)
            @test c_precal ≈ 0.01 rtol=1e-3
            @test c ≈ 0.01 rtol=1e-3
            @test (h_deconv.weights, positions) == peakfinder_result(h; threshold, σ=3.0)
        end

        h_uncal = RadiationSpectra.get_example_spectrum()
        c = calibrate_spectrum(h_uncal, lines; min_n_peaks=30, σ=2.0)[5]
        @test c ≈ 0.011259105696794384 rtol=1e-3
    end
end
