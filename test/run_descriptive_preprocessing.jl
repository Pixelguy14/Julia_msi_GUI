using CairoMakie
using Statistics
using Colors

function create_msi_parameter_plot()
    # Create the figure with subplots for each preprocessing step
    fig = Figure(size=(1600, 1200), fontsize=12)
    
    # Define the preprocessing steps and their parameters
    steps = [
        ("BaselineCorrection", ["iterations: 10", "method: SNIP", "window: 3.0"]),
        ("Calibration", ["fit_order: 2", "method: internal_standards", "ppm_tolerance: 20.0"]),
        ("Normalization", ["method: rms"]),
        ("PeakAlignment", ["max_shift_ppm: 50.0", "method: linear", "tolerance: 20.0 ppm"]),
        ("PeakBinningParams", ["max_bin_width_ppm: 60.0", "method: adaptive", "min_peak_per_bin: 3"]),
        ("PeakPicking", ["half_window: 5", "method: centroid", "snr_threshold: 15.0"]),
        ("PeakSelection", ["max_fwhm_ppm: 95.39", "min_fwhm_ppm: 19.08", "min_snr: 3.0"]),
        ("Smoothing", ["method: savitzky_golay", "order: 3", "window: 5"])
    ]
    
    # Create a grid of subplots
    g = fig[1, 1] = GridLayout()
    
    # Plot each preprocessing step
    for (idx, (step_name, params)) in enumerate(steps)
        row, col = fldmod1(idx, 2)
        ax = Axis(g[row, col], title=step_name, titlesize=14)
        
        # Generate simulated m/z values and intensities
        mz_range = range(100, 1000, length=500)
        
        if step_name == "BaselineCorrection"
            # Simulate spectrum with baseline
            true_signal = 0.5 .* exp.(-0.001 .* (mz_range .- 400).^2) .+
                          0.3 .* exp.(-0.002 .* (mz_range .- 600).^2)
            baseline = 0.1 .+ 0.05 .* sin.(mz_range ./ 50)
            noisy_signal = true_signal .+ baseline .+ 0.02 .* randn(length(mz_range))
            corrected = noisy_signal .- baseline
            
            lines!(ax, mz_range, noisy_signal, color=:blue, linewidth=2, label="Raw")
            lines!(ax, mz_range, baseline, color=:red, linewidth=2, linestyle=:dash, label="Baseline")
            lines!(ax, mz_range, corrected, color=:green, linewidth=2, label="Corrected")
            
        elseif step_name == "Calibration"
            # Simulate calibration shift
            reference_peaks = [200, 400, 600, 800]
            measured_peaks = reference_peaks .+ 2.0 .* randn(length(reference_peaks))
            
            scatter!(ax, reference_peaks, fill(0.5, length(reference_peaks)), 
                    color=:red, markersize=15, label="Reference")
            scatter!(ax, measured_peaks, fill(0.3, length(measured_peaks)), 
                    color=:blue, markersize=10, label="Measured")
            
            # Add calibration lines
            for i in 1:length(reference_peaks)
                lines!(ax, [measured_peaks[i], reference_peaks[i]], [0.3, 0.5], 
                      color=:black, linewidth=1, linestyle=:dash)
            end
            
        elseif step_name == "Normalization"
            # Simulate normalization effect
            spectra = [
                0.8 .* exp.(-0.001 .* (mz_range .- 300).^2) .+ 0.2 .* randn(length(mz_range)),
                1.2 .* exp.(-0.001 .* (mz_range .- 300).^2) .+ 0.2 .* randn(length(mz_range)),
                0.9 .* exp.(-0.001 .* (mz_range .- 300).^2) .+ 0.2 .* randn(length(mz_range))
            ]
            
            normalized_spectra = [spec ./ std(spec) for spec in spectra]
            
            for (i, spec) in enumerate(spectra)
                lines!(ax, mz_range, spec .+ i*0.3, color=RGBA(1, 0, 0, 0.6), linewidth=2, 
                      label=i==1 ? "Before Norm" : "")
            end
            for (i, spec) in enumerate(normalized_spectra)
                lines!(ax, mz_range, spec .+ i*0.3, color=RGBA(0, 0, 1, 0.6), linewidth=2, 
                      label=i==1 ? "After Norm" : "")
            end
            
        elseif step_name == "PeakAlignment"
            # Simulate peak alignment
            base_peaks = [300, 500, 700]
            shifts = [-15, 5, 10]
            
            for (i, shift) in enumerate(shifts)
                shifted_peaks = base_peaks .+ shift
                aligned_peaks = base_peaks
                
                scatter!(ax, shifted_peaks, fill(i, length(shifted_peaks)), 
                        color=:red, markersize=12, label=i==1 ? "Before Align" : "")
                scatter!(ax, aligned_peaks, fill(i+0.3, length(aligned_peaks)), 
                        color=:green, markersize=12, label=i==1 ? "After Align" : "")
                
                # Show alignment lines
                for j in 1:length(base_peaks)
                    lines!(ax, [shifted_peaks[j], aligned_peaks[j]], [i, i+0.3], 
                          color=:black, linewidth=1, linestyle=:dash)
                end
            end
            
        elseif step_name == "PeakBinningParams"
            # Simulate peak binning with centroids
            raw_peaks_mz = 100:25:900
            raw_peaks_intensity = rand(length(raw_peaks_mz))
            
            # Create binned peaks (wider bins)
            bin_centers = 150:60:850
            bin_intensities = [sum(raw_peaks_intensity[abs.(raw_peaks_mz .- center) .< 30]) 
                              for center in bin_centers] .* 0.8
            
            # Profile mode (continuous)
            profile_signal = zeros(length(mz_range))
            for (mz, int) in zip(raw_peaks_mz, raw_peaks_intensity)
                profile_signal .+= int .* exp.(-0.001 .* (mz_range .- mz).^2)
            end
            
            lines!(ax, mz_range, profile_signal, color=:blue, linewidth=2, label="Profile")
            barplot!(ax, bin_centers, bin_intensities, color=RGBA(1, 0, 0, 0.7), 
                    width=50, label="Binned Centroids")
            
        elseif step_name == "PeakPicking"
            # Simulate peak picking from profile to centroids
            profile_signal = 0.6 .* exp.(-0.0005 .* (mz_range .- 400).^2) .+
                             0.4 .* exp.(-0.0008 .* (mz_range .- 650).^2) .+
                             0.1 .* randn(length(mz_range))
            
            # Simulate picked peaks (centroids)
            peak_positions = [380, 405, 640, 660]
            peak_intensities = [0.5, 0.6, 0.35, 0.4]
            
            lines!(ax, mz_range, profile_signal, color=:blue, linewidth=2, label="Profile Spectrum")
            scatter!(ax, peak_positions, peak_intensities, color=:red, markersize=20, 
                    label="Picked Centroids", strokewidth=2)
            
        elseif step_name == "PeakSelection"
            # Simulate peak selection based on criteria
            all_peaks_mz = 200:50:800
            all_peaks_fwhm = rand(length(all_peaks_mz)) .* 100 .+ 10
            all_peaks_snr = rand(length(all_peaks_mz)) .* 10
            
            # Selection criteria
            selected = (all_peaks_fwhm .>= 19.08) .& (all_peaks_fwhm .<= 95.39) .& (all_peaks_snr .>= 3.0)
            
            scatter!(ax, all_peaks_mz[.!selected], all_peaks_fwhm[.!selected], 
                    color=:red, markersize=15, label="Rejected")
            scatter!(ax, all_peaks_mz[selected], all_peaks_fwhm[selected], 
                    color=:green, markersize=15, label="Selected")
            
            # Add selection criteria lines
            hlines!(ax, [19.08, 95.39], color=:black, linestyle=:dash, linewidth=2)
            text!(ax, 850, 50; text="FWHM bounds", color=:black, fontsize=10)
            
        elseif step_name == "Smoothing"
            # Simulate smoothing effect
            true_signal = 0.7 .* exp.(-0.001 .* (mz_range .- 450).^2) .+
                          0.5 .* exp.(-0.0008 .* (mz_range .- 650).^2)
            noisy_signal = true_signal .+ 0.1 .* randn(length(mz_range))
            
            # Simple smoothing simulation
            smoothed = similar(noisy_signal)
            window = 5
            for i in 1:length(noisy_signal)
                start_idx = max(1, i - window ÷ 2)
                end_idx = min(length(noisy_signal), i + window ÷ 2)
                smoothed[i] = mean(noisy_signal[start_idx:end_idx])
            end
            
            lines!(ax, mz_range, noisy_signal, color=:red, linewidth=1, label="Noisy")
            lines!(ax, mz_range, smoothed, color=:blue, linewidth=2, label="Smoothed")
            lines!(ax, mz_range, true_signal, color=:green, linewidth=1, linestyle=:dash, label="True")
        end
        
        # Add parameter text using a more reliable approach
        param_text = join(params, "\n")
        
        text!(ax, param_text, position=Point2f(0.05, 0.95), 
            space=:relative, align=(:left, :top), color=:black, 
            fontsize=10, font=:regular)
        
        # Add legend for selected plots
        if idx <= 4
            axislegend(ax, position=:rt, framevisible=true, backgroundcolor=RGBA(1,1,1,0.8))
        end
        
        # Customize axes
        ax.xlabel = "m/z"
        ax.ylabel = idx in [1,3,5,7] ? "Intensity" : ""
        ax.xgridvisible = false
        ax.ygridvisible = false
    end
    
    # Add overall title
    Label(fig[0, :], "MSI Preprocessing Pipeline Parameters and Simulations", 
          fontsize=18, font=:bold, padding=(0, 0, 10, 0))
    
    # Add explanation
    explanation = """
    Simulation of MSI preprocessing parameters showing:
    • Blue lines: Profile/continuous spectra
    • Red bars/points: Centroid data  
    • Dashed lines: Reference/true signals
    • Green: Processed/corrected data
    • Each subplot demonstrates key parameters for the preprocessing step
    """
    
    Label(fig[2, :], explanation, fontsize=12, tellwidth=false, padding=(10, 10, 10, 10))
    
    # Adjust layout
    colgap!(g, 20)
    rowgap!(g, 20)
    
    fig
end

# Create and display the plot
fig = create_msi_parameter_plot()

# Save the plot
save("msi_preprocessing_parameters.png", fig)
println("Plot saved as 'msi_preprocessing_parameters.png'")

# Display the plot (if in an interactive environment)
fig
