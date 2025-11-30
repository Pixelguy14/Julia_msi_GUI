using CairoMakie
using Statistics
using Colors

"""
Plots the effect of baseline correction, showing how different numbers of SNIP
iterations affect the estimated baseline.
"""
function plot_baseline_correction_details()
    fig = Figure(size=(1200, 700), fontsize=14)
    ax = Axis(fig[1, 1], title="Detailed View: Baseline Correction (SNIP)", xlabel="m/z", ylabel="Intensity")

    mz_range = range(200, 1000, length=1500)
    true_signal = 1.2 .* exp.(-0.001 .* (mz_range .- 450).^2) .+ 0.8 .* exp.(-0.002 .* (mz_range .- 750).^2)
    baseline_comp = 0.15 .+ 0.1 .* sin.(mz_range ./ 60) .+ 0.0001 .* mz_range
    noise = 0.04 .* randn(length(mz_range))
    raw_signal = true_signal .+ baseline_comp .+ noise

    # Simplified SNIP simulation for visualization
    function simple_snip(y, iterations)
        b = copy(y)
        for _ in 1:iterations
            for i in 2:length(b)-1
                b[i] = min(b[i], 0.5 * (b[i-1] + b[i+1]))
            end
        end
        return b
    end
    
    baseline_iter_20 = simple_snip(raw_signal, 20)
    baseline_iter_200 = simple_snip(raw_signal, 200)
    corrected_signal = raw_signal .- baseline_iter_200

    lines!(ax, mz_range, raw_signal, color=(:grey, 0.6), label="Raw Signal")
    lines!(ax, mz_range, corrected_signal, color=:green, linewidth=2.5, label="Corrected Signal")
    l1 = lines!(ax, mz_range, baseline_iter_20, color=(:red, 0.7), linestyle=:dash, linewidth=2, label="Baseline (iterations: 20)")
    l2 = lines!(ax, mz_range, baseline_iter_200, color=:red, linewidth=2.5, label="Baseline (iterations: 200)")

    text!(ax, "Parameter: `iterations`\nMore iterations result in a more aggressive baseline that follows the signal floor more closely.", 
          position=Point2f(0.05, 0.95), space=:relative, align=(:left, :top), fontsize=12)
    axislegend(ax, position=:rt)
    save("descriptive_plot_baseline_correction.png", fig)
    println("Saved: descriptive_plot_baseline_correction.png")
end

"""
Plots the effect of smoothing, comparing different window sizes.
"""
function plot_smoothing_details()
    fig = Figure(size=(1200, 700), fontsize=14)
    ax = Axis(fig[1, 1], title="Detailed View: Smoothing (Savitzky-Golay)", xlabel="m/z", ylabel="Intensity")

    mz_range = range(400, 700, length=1000)
    true_signal = 0.8 .* exp.(-0.002 .* (mz_range .- 500).^2) .+ 0.6 .* exp.(-0.001 .* (mz_range .- 600).^2)
    noisy_signal = true_signal .+ 0.1 .* randn(length(mz_range))

    function simple_moving_average(y, window)
        smoothed = similar(y)
        for i in 1:length(y)
            start_idx = max(1, i - window ÷ 2)
            end_idx = min(length(y), i + window ÷ 2)
            smoothed[i] = mean(y[start_idx:end_idx])
        end
        return smoothed
    end
    
    smoothed_small_window = simple_moving_average(noisy_signal, 5)
    smoothed_large_window = simple_moving_average(noisy_signal, 21)

    lines!(ax, mz_range, noisy_signal, color=(:red, 0.4), label="Noisy Signal")
    lines!(ax, mz_range, true_signal, color=:black, linestyle=:dash, linewidth=2, label="True Signal")
    lines!(ax, mz_range, smoothed_small_window, color=:blue, linewidth=2, label="Smoothed (window: 5)")
    lines!(ax, mz_range, smoothed_large_window, color=:purple, linewidth=2, label="Smoothed (window: 21)")
    
    text!(ax, "Parameter: `window`\nA larger window increases smoothing but may broaden peaks.", 
          position=Point2f(0.05, 0.95), space=:relative, align=(:left, :top), fontsize=12)
    axislegend(ax, position=:rt)
    save("descriptive_plot_smoothing.png", fig)
    println("Saved: descriptive_plot_smoothing.png")
end

"""
Visualizes the peak picking process for profile-mode data, illustrating the
effects of SNR threshold and peak prominence.
"""
function plot_peak_picking_profile_details()
    fig = Figure(size=(1200, 700), fontsize=14)
    ax = Axis(fig[1, 1], title="Detailed View: Peak Picking (Profile Mode)", xlabel="m/z", ylabel="Intensity")

    mz = 1:200
    base_signal = 10 .* exp.(-((mz .- 50).^2) ./ (2*3^2)) .+ 7 .* exp.(-((mz .- 120).^2) ./ (2*5^2)) .+ 2
    small_peak_signal = zeros(200)
    for i in 80:90
        small_peak_signal[i] = 3 * exp(-((i - 85)^2) / 2.0)
    end
    noise = 0.5 .* randn(200)
    intensity = base_signal .+ small_peak_signal .+ noise
    
    noise_level = median(abs.(intensity .- median(intensity))) * 1.4826 # MAD
    snr_threshold_val = 3.0
    intensity_threshold = noise_level * snr_threshold_val

    picked_peaks_mz = [50, 85, 120]
    picked_peaks_intensity = intensity[picked_peaks_mz]
    
    hlines!(ax, [noise_level], color=:gray, linestyle=:dot, label="Est. Noise Level")
    hlines!(ax, [intensity_threshold], color=:orange, linestyle=:dash, label="SNR Threshold (snr_threshold = 3.0)")
    lines!(ax, mz, intensity, color=:blue, label="Profile Spectrum")
    scatter!(ax, picked_peaks_mz, picked_peaks_intensity, color=:green, markersize=15, strokewidth=2, label="Peaks passing SNR")
    scatter!(ax, [25], [intensity[25]], color=:red, marker=:x, markersize=15, label="Local max below SNR")

    # Illustrate prominence
    arrows!(ax, [120, 120], [intensity[135], intensity[120]], [0, 0], [intensity[120]-intensity[135], 0], color=:purple)
    text!(ax, 125, (intensity[120]+intensity[135])/2, text="Prominence", color=:purple)

    axislegend(ax, position=:rt)
    save("descriptive_plot_peakpicking_profile.png", fig)
    println("Saved: descriptive_plot_peakpicking_profile.png")
end

"""
Visualizes peak picking (filtering) for centroid-mode data based on an SNR
(intensity) threshold.
"""
function plot_peak_picking_centroid_details()
    fig = Figure(size=(1200, 700), fontsize=14)
    ax = Axis(fig[1, 1], title="Detailed View: Peak Picking (Centroid Mode)", xlabel="m/z", ylabel="Intensity")

    mz = [100, 150, 200, 250, 300, 350, 400]
    intensity = [10, 5, 25, 8, 3, 18, 12]
    snr_threshold_val = 10.0

    stem!(ax, mz, intensity, color=:gray, label="Input Centroids")
    
    selected_mask = intensity .>= snr_threshold_val
    stem!(ax, mz[selected_mask], intensity[selected_mask], color=:green, trunkwidth=3, label="Selected Peaks (intensity >= 10)")
    stem!(ax, mz[.!selected_mask], intensity[.!selected_mask], color=:red, trunkwidth=3, label="Rejected Peaks (intensity < 10)")
    
    hlines!(ax, [snr_threshold_val], color=:orange, linestyle=:dash, label="Intensity Threshold (snr_threshold)")

    text!(ax, "Parameter: `snr_threshold`\nIn centroid mode, this acts as a direct intensity filter.", 
          position=Point2f(0.05, 0.95), space=:relative, align=(:left, :top), fontsize=12)
    axislegend(ax, position=:rt)
    save("descriptive_plot_peakpicking_centroid.png", fig)
    println("Saved: descriptive_plot_peakpicking_centroid.png")
end

"""
Plots the effect of different normalization methods on a set of spectra.
"""
function plot_normalization_details()
    fig = Figure(size=(1200, 700), fontsize=14)
    
    mz = range(300, 400, length=500)
    spec1 = 1.5 .* exp.(-((mz .- 350).^2) ./ 50)
    spec2 = 0.8 .* exp.(-((mz .- 350).^2) ./ 50)
    
    ax1 = Axis(fig[1, 1], title="Before Normalization", ylabel="Absolute Intensity")
    lines!(ax1, mz, spec1, label="Spectrum A (High TIC)")
    lines!(ax1, mz, spec2, label="Spectrum B (Low TIC)")
    axislegend(ax1)

    ax2 = Axis(fig[1, 2], title="After Normalization (TIC)", ylabel="Relative Intensity")
    lines!(ax2, mz, spec1 ./ sum(spec1), label="Spectrum A (Normalized)")
    lines!(ax2, mz, spec2 ./ sum(spec2), label="Spectrum B (Normalized)")

    text!(ax2, "Effect: Spectra are scaled to have the same total area, making their intensities comparable.",
          position=Point2f(0.05, 0.95), space=:relative, align=(:left, :top), fontsize=12, justification=:left)
    
    save("descriptive_plot_normalization.png", fig)
    println("Saved: descriptive_plot_normalization.png")
end

"""
Illustrates mass calibration, showing how a calibration curve corrects measured
m/z values based on reference peaks.
"""
function plot_calibration_details()
    fig = Figure(size=(1200, 700), fontsize=14)
    ax = Axis(fig[1, 1], title="Detailed View: Calibration", xlabel="Measured m/z", ylabel="m/z Error (Measured - Reference)")

    ref_mz = [200, 400, 600, 800, 1000]
    measured_mz = ref_mz .+ [0.1, 0.15, 0.2, 0.25, 0.3] .+ 0.02 .* randn(5)
    errors = measured_mz .- ref_mz

    # Fit a linear model to the error
    A = [ones(5) measured_mz]
    coeffs = A \ errors
    correction_func(m) = m - (coeffs[1] .+ coeffs[2] .* m)
    
    fit_line = coeffs[1] .+ coeffs[2] .* measured_mz
    
    scatter!(ax, measured_mz, errors, color=:red, markersize=15, label="Measured Error")
    lines!(ax, measured_mz, fit_line, color=:blue, label="Calibration Curve (fit_order=1)")

    text!(ax, "Parameter: `fit_order`\nA curve is fit to the error of known reference peaks.\nThis curve is then used to correct all m/z values.", 
          position=Point2f(0.05, 0.95), space=:relative, align=(:left, :top), fontsize=12)
    axislegend(ax, position=:rb)
    save("descriptive_plot_calibration.png", fig)
    println("Saved: descriptive_plot_calibration.png")
end

"""
Visualizes peak alignment by showing multiple spectra with misaligned peaks
before and after the alignment process.
"""
function plot_alignment_details()
    fig = Figure(size=(1600, 600), fontsize=14)
    ax1 = Axis(fig[1, 1], title="Before Alignment", xlabel="m/z", yticklabelsvisible=false, ygridvisible=false)
    ax2 = Axis(fig[1, 2], title="After Alignment", xlabel="m/z", yticklabelsvisible=false, ygridvisible=false)

    mz = range(490, 510, length=1000)
    shifts = [-0.5, 0.0, 0.8]
    colors = [:blue, :green, :purple]
    
    for (i, shift) in enumerate(shifts)
        peak_center = 500 + shift
        spectrum = exp.(-((mz .- peak_center).^2) ./ 0.1)
        lines!(ax1, mz, spectrum .+ i, color=colors[i])
        vlines!(ax1, [peak_center], color=(colors[i], 0.5), linestyle=:dash)

        # After alignment, all peaks are at 500
        aligned_spectrum = exp.(-((mz .- 500).^2) ./ 0.1)
        lines!(ax2, mz, aligned_spectrum .+ i, color=colors[i])
    end
    vlines!(ax2, [500], color=:red, linestyle=:dash, label="Reference m/z")
    axislegend(ax2)

    save("descriptive_plot_alignment.png", fig)
    println("Saved: descriptive_plot_alignment.png")
end

"""
Illustrates peak selection by filtering a population of peaks based on
FWHM and SNR criteria.
"""
function plot_peak_selection_details()
    fig = Figure(size=(1200, 700), fontsize=14)
    ax = Axis(fig[1, 1], title="Detailed View: Peak Selection", xlabel="FWHM (ppm)", ylabel="Signal-to-Noise Ratio (SNR)")

    n_peaks = 100
    fwhm = rand(n_peaks) .* 150
    snr = rand(n_peaks) .* 20
    
    min_fwhm_ppm = 20.0
    max_fwhm_ppm = 100.0
    min_snr = 5.0

    selected_mask = (fwhm .>= min_fwhm_ppm) .& (fwhm .<= max_fwhm_ppm) .& (snr .>= min_snr)

    scatter!(ax, fwhm[.!selected_mask], snr[.!selected_mask], color=(:red, 0.5), label="Rejected Peaks")
    scatter!(ax, fwhm[selected_mask], snr[selected_mask], color=:green, label="Selected Peaks")

    vlines!(ax, [min_fwhm_ppm, max_fwhm_ppm], color=:blue, linestyle=:dash, label="FWHM bounds")
    hlines!(ax, [min_snr], color=:orange, linestyle=:dash, label="SNR bound")
    
    poly!(ax, BBox(min_fwhm_ppm, max_fwhm_ppm, min_snr, 22), color=(:green, 0.1))
    text!(ax, "Selection Region", position=(60, 12), color=:green, fontsize=14)

    axislegend(ax)
    save("descriptive_plot_peak_selection.png", fig)
    println("Saved: descriptive_plot_peak_selection.png")
end

"""
Visualizes the adaptive peak binning process, showing how peaks from different
spectra are grouped into a common bin based on a PPM tolerance.
"""
function plot_peak_binning_details()
    fig = Figure(size=(1200, 700), fontsize=14)
    ax = Axis(fig[1, 1], title="Detailed View: Adaptive Peak Binning", xlabel="m/z", yticklabelsvisible=false)

    ref_mz = 500.0
    tolerance_ppm = 50.0
    tol_mz = ref_mz * tolerance_ppm / 1e6
    
    bin_start = ref_mz - tol_mz/2
    bin_end = ref_mz + tol_mz/2
    
    peaks_mz = [ref_mz - 0.01, ref_mz + 0.005, ref_mz + 0.02, ref_mz - 0.015]
    peak_intensities = [0.8, 1.0, 0.9, 0.7]
    peak_colors = [:blue, :green, :purple, :orange]
    
    vspan!(ax, bin_start, bin_end, color=(:gray, 0.2), label="Bin (tolerance: 50 ppm)")
    stem!(ax, peaks_mz, peak_intensities, color=peak_colors, markersize=15)
    
    # Show bin center
    bin_center = mean(peaks_mz)
    vlines!(ax, [bin_center], color=:red, linestyle=:dash, label="Calculated Bin Center")

    text!(ax, "Parameter: `tolerance`\nPeaks from different spectra within the tolerance window are grouped into a single feature.",
          position=Point2f(0.05, 0.95), space=:relative, align=(:left, :top), fontsize=12)
    axislegend(ax, position=:rt)
    xlims!(ax, ref_mz - tol_mz*2, ref_mz + tol_mz*2)
    save("descriptive_plot_peak_binning.png", fig)
    println("Saved: descriptive_plot_peak_binning.png")
end


"""
Main function to generate and save all descriptive plots.
"""
function create_and_save_all_plots()
    println("Generating detailed descriptive plots for preprocessing steps...")
    
    plot_baseline_correction_details()
    plot_smoothing_details()
    plot_peak_picking_profile_details()
    plot_peak_picking_centroid_details()
    plot_normalization_details()
    plot_calibration_details()
    plot_alignment_details()
    plot_peak_selection_details()
    plot_peak_binning_details()

    println("\nAll descriptive plots have been saved in the current directory.")
end

# Execute the plot generation
if abspath(PROGRAM_FILE) == @__FILE__
    create_and_save_all_plots()
end

