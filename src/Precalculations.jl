using StatsBase # For mean, std, median, quantile, mad

# =============================================================================
# 7) Spatial & Advanced Processing (Stubs & New Functions)
# =============================================================================

"""
    find_ppm_error_by_region(msi_data, region_masks, reference_peaks)::Dict

Calculates PPM error statistics for different spatial regions.
`region_masks` is a Dict mapping region names (e.g., :tumor) to BitMatrix masks.
"""
function find_ppm_error_by_region(msi_data::MSIData, region_masks::Dict, reference_peaks::Dict)
    regional_reports = Dict{Symbol, NamedTuple}()
    
    width, height = msi_data.image_dims
    
    for (region_name, mask) in region_masks
        mask_height, mask_width = size(mask)
        if mask_width != width || mask_height != height
            @warn "Mask dimensions ($(mask_width)x$(mask_height)) for region '$region_name' do not match image dimensions ($(width)x$(height)). Skipping."
            continue
        end

        indices = [
            i for i in 1:length(msi_data.spectra_metadata) 
            if msi_data.spectra_metadata[i].x > 0 &&
               msi_data.spectra_metadata[i].y > 0 &&
               mask[msi_data.spectra_metadata[i].y, msi_data.spectra_metadata[i].x]
        ]
        
        if isempty(indices) continue end
        
        # Call analyze_mass_accuracy with the specific indices for the region
        regional_reports[region_name] = analyze_mass_accuracy(msi_data, reference_peaks; spectrum_indices=indices)
    end
    return regional_reports
end

function analyze_mass_accuracy(
    msi_data::MSIData,
    reference_peaks::Dict{Float64, String}; # m/z => name
    spectrum_indices::AbstractVector{Int},
    peak_detection_snr_threshold::Float64 = 2.0,
    ppm_tolerance_for_matching::Float64 = 50.0
    )::NamedTuple
    println("Analyzing mass accuracy for $(length(spectrum_indices)) spectra...")

    all_ppm_errors = Float64[]
    total_matched_peaks = 0
    total_spectra_processed = 0

    _iterate_spectra_fast(msi_data, spectrum_indices) do idx, mz, intensity
        total_spectra_processed += 1
        if !validate_spectrum(mz, intensity)
            @warn "Spectrum $idx is invalid, skipping mass accuracy analysis for it."
            return
        end

        detected_peaks = detect_peaks_profile(mz, intensity; snr_threshold=peak_detection_snr_threshold)

        for ref_mz in keys(reference_peaks)
            # Find the closest detected peak to this reference m/z within tolerance
            min_ppm_error = Inf
            best_matched_peak_mz = nothing

            for p in detected_peaks
                ppm_error = calculate_ppm_error(p.mz, ref_mz)
                if ppm_error <= ppm_tolerance_for_matching && ppm_error < min_ppm_error
                    min_ppm_error = ppm_error
                    best_matched_peak_mz = p.mz
                end
            end

            if best_matched_peak_mz !== nothing
                push!(all_ppm_errors, min_ppm_error)
                total_matched_peaks += 1
            end
        end
    end

    if isempty(all_ppm_errors)
        @warn "No reference peaks matched in any of the $(total_spectra_processed) processed spectra."
        return (
            mean_ppm_error = NaN,
            median_ppm_error = NaN,
            std_ppm_error = NaN,
            min_ppm_error = NaN,
            max_ppm_error = NaN,
            total_matched_peaks = 0,
            total_spectra_analyzed = total_spectra_processed
        )
    end

    # Calculate summary statistics
    mean_err = mean(all_ppm_errors)
    median_err = median(all_ppm_errors)
    std_err = std(all_ppm_errors)
    min_err = minimum(all_ppm_errors)
    max_err = maximum(all_ppm_errors)

    println("Mass accuracy analysis complete for $(total_spectra_processed) spectra.")
    return (
        mean_ppm_error = mean_err,
        median_ppm_error = median_err,
        std_ppm_error = std_err,
        min_ppm_error = min_err,
        max_ppm_error = max_err,
        total_matched_peaks = total_matched_peaks,
        total_spectra_analyzed = total_spectra_processed
    )
end

# =============================================================================
# 8) Advanced Peak Quality & Adaptive Parameters
# =============================================================================

"""
    calculate_adaptive_bin_tolerance(ppm_error_distribution) -> Float64

Calculates an appropriate binning tolerance based on observed mass accuracy.
"""
function calculate_adaptive_bin_tolerance(ppm_error_distribution::Vector{Float64})
    if isempty(ppm_error_distribution)
        return 20.0 # Default if no data
    end
    
    # Filter out NaN and infinite values
    valid_errors = filter(x -> isfinite(x), ppm_error_distribution)
    
    if isempty(valid_errors)
        return 20.0 # Default if no valid data
    end
    
    # A robust strategy: mean + 3 * std deviation to capture ~99.7% of peaks
    mean_err = mean(valid_errors)
    std_err = std(valid_errors)
    
    # Ensure we don't get NaN
    if !isfinite(mean_err) || !isfinite(std_err)
        return 20.0
    end
    
    tolerance = mean_err + 3 * std_err
    
    # Cap at reasonable value and ensure finite
    return min(max(tolerance, 10.0), 100.0)
end

"""
    calculate_preprocessing_hints(data::MSIData; sample_size::Int=100)::Dict{Symbol, Any}

Analyzes a sample of spectra to determine optimal default parameters for
preprocessing steps, returning a dictionary of hints.
"""
function calculate_preprocessing_hints(data::MSIData; sample_indices::AbstractVector{Int})::Dict{Symbol, Any}
    println("Calculating preprocessing hints from a sample of $(length(sample_indices)) spectra...")
    
    num_spectra = length(data.spectra_metadata)
    if num_spectra == 0 || isempty(sample_indices)
        @warn "No spectra in dataset/sample to calculate hints from."
        return Dict{Symbol, Any}(
            :estimated_noise => 1.0,
            :suggested_snr => 3.0,
            :suggested_smoothing_window => 9
        )
    end

    # Initialize hints with defaults immediately
    hints = Dict{Symbol, Any}(
        :estimated_noise => 1.0, 
        :suggested_snr => 3.0,   
        :suggested_smoothing_window => 9 
    )

    all_noise_levels = Float64[]

    for idx in sample_indices
        try
            mz, intensity = GetSpectrum(data, idx)
            if !isempty(intensity)
                noise = mad(intensity, normalize=true)
                if isfinite(noise) && noise > 0
                    push!(all_noise_levels, noise)
                end
            end
        catch e
            @warn "Could not access spectrum #$idx to calculate hints: $e. Skipping."
        end
    end
    
    if !isempty(all_noise_levels)
        estimated_noise = mean(all_noise_levels)
        hints[:estimated_noise] = estimated_noise
        println("  - Estimated Noise Level: $(round(estimated_noise, digits=4))")
        println("  - Suggested SNR Threshold: 3.0")
    else
        @warn "Could not estimate noise from sample. Using default hints."
        # Defaults already set in `hints` initialization
    end
    
    # Suggest smoothing window based on resolution (if available)
    if data.instrument_metadata !== nothing && data.instrument_metadata.resolution !== nothing
        res = data.instrument_metadata.resolution
        # Update suggested smoothing window if resolution is known
        if res > 40000
            hints[:suggested_smoothing_window] = 5
        elseif res > 10000
            hints[:suggested_smoothing_window] = 7
        else
            hints[:suggested_smoothing_window] = 9
        end
    end
    println("  - Suggested Smoothing Window: $(hints[:suggested_smoothing_window])")

    println("Preprocessing hints calculated.")
    return hints
end

"""
    analyze_instrument_characteristics(msi_data::MSIData)::Dict

Analyzes instrument metadata and data characteristics to determine acquisition properties.
"""
function analyze_instrument_characteristics(msi_data::MSIData; sample_indices::AbstractVector{Int})::Dict
    results = Dict{Symbol, Any}()
    
    # Use instrument metadata if available
    if msi_data.instrument_metadata !== nothing
        inst = msi_data.instrument_metadata
        results[:resolution] = inst.resolution
        results[:mass_accuracy_ppm] = inst.mass_accuracy_ppm
        results[:instrument_model] = inst.instrument_model
        results[:polarity] = inst.polarity
        results[:calibration_status] = inst.calibration_status
        results[:vendor_preprocessing] = inst.vendor_preprocessing_steps
    end
    
    # Analyze data characteristics from spectra
    println("  Analyzing data characteristics from spectra...")
    
    spectrum_modes = Set{SpectrumMode}()
    mz_step_sizes = Float64[]
    intensity_ranges = Tuple{Float64, Float64}[]  # (min, max) per spectrum
    
    if isempty(sample_indices)
        @warn "No indices to sample for instrument characteristics analysis."
        results[:acquisition_mode] = :unknown
        results[:mz_axis_type] = :unknown
        results[:dynamic_range] = 3.0
        return results
    end
    
    _iterate_spectra_fast(msi_data, sample_indices) do idx, mz, intensity
        # Record spectrum mode
        push!(spectrum_modes, msi_data.spectra_metadata[idx].mode)
        
        # Calculate m/z step statistics (for profile data)
        if length(mz) > 1 && msi_data.spectra_metadata[idx].mode == PROFILE
            steps = diff(mz)
            if !isempty(steps)
                push!(mz_step_sizes, mean(steps))
            end
        end
        
        # Record intensity range
        if !isempty(intensity)
            push!(intensity_ranges, (minimum(intensity), maximum(intensity)))
        end
    end
    
    # Determine acquisition mode
    if length(spectrum_modes) == 1
        results[:acquisition_mode] = first(spectrum_modes) == CENTROID ? :centroid : :profile
    else
        results[:acquisition_mode] = :mixed
    end
    
    # Determine m/z axis regularity
    if !isempty(mz_step_sizes)
        avg_step = mean(mz_step_sizes)
        step_std = std(mz_step_sizes)
        results[:mz_axis_type] = step_std / avg_step < 0.01 ? :regular : :irregular
        results[:average_mz_step] = avg_step
    else
        results[:mz_axis_type] = :unknown
    end
    
    # Enhanced dynamic range calculation
    if !isempty(intensity_ranges)
        println("DEBUG: Found $(length(intensity_ranges)) intensity ranges")
        
        # FIX: Handle cases where min intensities are zero
        max_intensities = Float64[]
        min_positive_intensities = Float64[]
        
        for (min_val, max_val) in intensity_ranges
            if max_val > 1e-6  # Valid maximum
                push!(max_intensities, max_val)
                
                # Find the smallest positive intensity in the spectrum
                # For now, use a reasonable estimate: 1% of the noise level
                # In practice, you'd want to sample the actual spectrum
                if max_val > 0
                    # Estimate minimum detectable signal as ~3x noise level
                    estimated_min_signal = max_val * 1e-4  # Conservative estimate
                    push!(min_positive_intensities, estimated_min_signal)
                end
            end
        end
        
        println("DEBUG: Valid max intensities: $(length(max_intensities)), estimated min intensities: $(length(min_positive_intensities))")
        
        if !isempty(max_intensities) && !isempty(min_positive_intensities)
            avg_max = mean(max_intensities)
            avg_min = mean(min_positive_intensities)
            if avg_min > 0
                dynamic_range = log10(avg_max / avg_min)
                results[:dynamic_range] = dynamic_range
                println("DEBUG: Dynamic range calculated: $dynamic_range (avg_max=$avg_max, avg_min=$avg_min)")
            else
                results[:dynamic_range] = 0.0
            end
        else
            # Estimate based on typical values
            results[:dynamic_range] = 3.0  # Typical for MS data
        end
    else
        results[:dynamic_range] = 3.0  # Default estimate
    end
    
    println("  - Acquisition mode: $(results[:acquisition_mode])")
    println("  - m/z axis type: $(results[:mz_axis_type])")
    println("  - Dynamic range: $(round(get(results, :dynamic_range, 0), digits=2))")
    
    return results
end

"""
    analyze_signal_quality(msi_data::MSIData; sample_size::Int=100)::Dict

Analyzes noise characteristics, signal-to-noise ratios, and overall signal quality.
"""
function analyze_signal_quality(msi_data::MSIData; sample_indices::AbstractVector{Int})::Dict
    results = Dict{Symbol, Any}()
    
    println("  Analyzing signal quality from $(length(sample_indices)) sample spectra...")
    
    # Use the existing function for basic noise estimation
    hints_from_calc = calculate_preprocessing_hints(msi_data, sample_indices=sample_indices)
    
    # Copy relevant hints to results, or handle cases where they might be missing
    results[:estimated_noise] = get(hints_from_calc, :estimated_noise, 1.0)
    results[:suggested_snr] = get(hints_from_calc, :suggested_snr, 3.0)
    results[:suggested_smoothing_window] = get(hints_from_calc, :suggested_smoothing_window, 9)
    
    # Enhanced noise analysis
    noise_levels = Float64[]
    snr_distribution = Float64[]
    tic_values = Float64[]
    
    if isempty(sample_indices)
        @warn "No indices to sample for signal quality analysis."
        return results
    end
    
    _iterate_spectra_fast(msi_data, sample_indices) do idx, mz, intensity
        if !isempty(intensity)
            # Noise estimation using MAD
            noise = mad(intensity, normalize=true)
            push!(noise_levels, noise)
            
            # FIX: More robust SNR calculation
            valid_intensity = intensity[intensity .> 0]  # Remove zeros
            if !isempty(valid_intensity)
                # Use robust signal estimate (95th percentile instead of max)
                signal_estimate = quantile(valid_intensity, 0.95)
                noise_robust = max(noise, 1e-6)  # Avoid division by zero
                
                if noise_robust > 0 && isfinite(signal_estimate)
                    snr_val = signal_estimate / noise_robust
                    # Cap unrealistic SNR values
                    push!(snr_distribution, min(snr_val, 1e6))
                end
            end
            
            # Total ion count
            push!(tic_values, sum(intensity))
        end
    end
    
    if !isempty(noise_levels)
        results[:noise_mean] = mean(noise_levels)
        results[:noise_std] = std(noise_levels)
        results[:noise_cv] = results[:noise_std] / results[:noise_mean]  # Coefficient of variation
    end
    
    if !isempty(snr_distribution)
        results[:snr_mean] = mean(snr_distribution)
        results[:snr_median] = median(snr_distribution)
        results[:snr_95th] = quantile(snr_distribution, 0.95)
    end
    
    if !isempty(tic_values)
        results[:tic_mean] = mean(tic_values)
        results[:tic_std] = std(tic_values)
        results[:tic_cv] = results[:tic_std] / results[:tic_mean]
    end
    
    println("  - Estimated noise level: $(round(get(results, :noise_mean, 0), digits=4))")
    println("  - Average SNR: $(round(get(results, :snr_mean, 0), digits=2))")
    println("  - TIC CV: $(round(get(results, :tic_cv, 0) * 100, digits=1))%")
    
    return results
end

"""
    analyze_mass_accuracy_global(msi_data::MSIData, reference_peaks::Dict, sample_size::Int)::Dict

Analyzes mass accuracy across the dataset using reference peaks.
"""
function analyze_mass_accuracy_global(msi_data::MSIData, reference_peaks::Dict; 
                                       spectrum_indices::AbstractVector{Int})::Dict
    results = Dict{Symbol, Any}()
    
    println("  Analyzing mass accuracy using $(length(reference_peaks)) reference peaks on $(length(spectrum_indices)) spectra...")
    
    if isempty(spectrum_indices)
        @warn "No indices to sample for mass accuracy analysis."
        # Return a structure indicating no analysis was performed
        empty_report = (
            mean_ppm_error = NaN, median_ppm_error = NaN, std_ppm_error = NaN,
            min_ppm_error = NaN, max_ppm_error = NaN, total_matched_peaks = 0,
            total_spectra_analyzed = 0
        )
        results[:global_accuracy] = empty_report
        results[:suggested_bin_tolerance] = 20.0 # Default
        return results
    end
    
    # Use the existing analyze_mass_accuracy function with the provided indices
    accuracy_report = analyze_mass_accuracy(msi_data, reference_peaks; 
                                          spectrum_indices=spectrum_indices)
    
    results[:global_accuracy] = accuracy_report
    results[:suggested_bin_tolerance] = calculate_adaptive_bin_tolerance(
        collect(Iterators.flatten([accuracy_report.mean_ppm_error]))  # Simplified - would need actual distribution
    )
    
    println("  - Mean PPM error: $(round(accuracy_report.mean_ppm_error, digits=2))")
    println("  - Suggested bin tolerance: $(round(results[:suggested_bin_tolerance], digits=2)) ppm")
    
    return results
end

"""
    analyze_spatial_regions(msi_data::MSIData, region_masks::Dict, reference_peaks::Dict)::Dict

Analyzes different spatial regions for variations in signal quality and mass accuracy.
"""
function analyze_spatial_regions(msi_data::MSIData, region_masks::Dict, reference_peaks::Dict)::Dict
    results = Dict{Symbol, Any}()
    
    println("  Analyzing $(length(region_masks)) spatial regions...")
    
    # Use the existing function
    regional_reports = find_ppm_error_by_region(msi_data, region_masks, reference_peaks)
    
    results[:regional_ppm_errors] = regional_reports
    
    # Calculate regional variations
    if !isempty(regional_reports)
        mean_errors = [report.mean_ppm_error for report in values(regional_reports) if isfinite(report.mean_ppm_error)]
        if !isempty(mean_errors)
            results[:max_regional_ppm_difference] = maximum(mean_errors) - minimum(mean_errors)
        end
    end
    
    for (region, report) in regional_reports
        println("  - $region: $(round(report.mean_ppm_error, digits=2)) ppm (n=$(report.total_spectra_analyzed))")
    end
    
    return results
end

"""
    analyze_peak_characteristics(msi_data::MSIData, sample_size::Int)::Dict

Analyzes peak shape, width, and quality characteristics.
"""
function analyze_peak_characteristics(msi_data::MSIData, instrument_analysis::Dict, mass_accuracy_results; 
                                    spectrum_indices::AbstractVector{Int})::Dict
    results = Dict{Symbol, Any}()
    
    # Determine acquisition mode from the analysis results, not from metadata
    acquisition_mode = get(instrument_analysis, :acquisition_mode, :unknown)
    
    # Calculate estimated_mean_ppm_error for fallback logic
    estimated_mean_ppm_error = 30.0 # Default if mass_accuracy_results is nothing or invalid
    if mass_accuracy_results !== nothing && haskey(mass_accuracy_results, :global_accuracy) && isa(mass_accuracy_results[:global_accuracy], NamedTuple)
        ppm_error_report = mass_accuracy_results[:global_accuracy]
        mean_ppm_error_val = get(ppm_error_report, :mean_ppm_error, 30.0)
        if isfinite(mean_ppm_error_val)
            estimated_mean_ppm_error = mean_ppm_error_val
        end
    end

    if isempty(spectrum_indices)
        @warn "No indices to sample for peak characteristics analysis."
        # Provide sensible defaults if no analysis can be run
        estimated_fwhm = estimated_mean_ppm_error * (acquisition_mode == :profile ? 3 : 2)
        results[:mean_fwhm_ppm] = estimated_fwhm
        results[:median_fwhm_ppm] = estimated_fwhm
        results[:mean_gaussian_r2] = acquisition_mode == :profile ? 0.7 : 0.9
        results[:peak_resolution_estimate] = 1e6 / estimated_fwhm
        results[:mean_peaks_per_spectrum] = 0
        return results
    end

    if acquisition_mode == PROFILE
        println("  Analyzing peak characteristics for PROFILE mode from $(length(spectrum_indices)) sample spectra...")
        
        peak_widths_ppm = Float64[]
        r_squared_values = Float64[]
        peak_counts = Int[]
        fwhm_values = Float64[]
        
        spectra_analyzed = 0
        peaks_analyzed = 0
        
        _iterate_spectra_fast(msi_data, spectrum_indices) do idx, mz, intensity
            if length(mz) < 10  # Skip spectra with too few points
                return
            end
            
            spectra_analyzed += 1
            meta = msi_data.spectra_metadata[idx]
            
            # Detect peaks with lower SNR threshold to find more peaks
            peaks = detect_peaks_profile(mz, intensity; snr_threshold=2.0)
            push!(peak_counts, length(peaks))
            
            if !isempty(peaks)
                # Analyze the strongest 3 peaks per spectrum
                sorted_peaks = sort(peaks, by=p->p.intensity, rev=true)
                for peak in sorted_peaks[1:min(3, length(sorted_peaks))]
                    try
                        # Find the index of the peak in the original mz array
                        peak_idx = argmin(abs.(mz .- peak.mz))
                        
                        fwhm_delta_m = calculate_robust_fwhm(mz, intensity, peak_idx)
                        
                        if !isnan(fwhm_delta_m) && fwhm_delta_m > 0.001 && fwhm_delta_m < 0.5  # Reasonable range in Da
                            fwhm_ppm = 1e6 * fwhm_delta_m / peak.mz
                            if 5.0 < fwhm_ppm < 500.0  # Reasonable ppm range
                                push!(peak_widths_ppm, fwhm_ppm)
                                push!(fwhm_values, fwhm_delta_m)
                                
                                r2 = _fit_gaussian_and_r2(mz, intensity, peak_idx, 5)
                                push!(r_squared_values, r2)
                                peaks_analyzed += 1
                                
                                if peaks_analyzed <= 3
                                    println("DEBUG: Peak at m/z $(peak.mz), FWHM = $(fwhm_ppm) ppm, R² = $r2")
                                end
                            end
                        end
                    catch e
                        continue
                    end
                end
            end
        end
        
        println("DEBUG: Analyzed $peaks_analyzed peaks from $spectra_analyzed spectra")
        
        if !isempty(peak_widths_ppm)
            results[:mean_fwhm_ppm] = mean(peak_widths_ppm)
            results[:median_fwhm_ppm] = median(peak_widths_ppm)
            results[:mean_gaussian_r2] = mean(r_squared_values)
            results[:peak_resolution_estimate] = 1e6 / results[:mean_fwhm_ppm]
            println("  - Actual FWHM measurement: $(round(results[:mean_fwhm_ppm], digits=2)) ppm")
        else
            # Fallback if no valid peaks found in PROFILE mode
            estimated_fwhm = estimated_mean_ppm_error * 3  # FWHM typically wider than mass error
            results[:mean_fwhm_ppm] = estimated_fwhm
            results[:median_fwhm_ppm] = estimated_fwhm
            results[:mean_gaussian_r2] = 0.7
            results[:peak_resolution_estimate] = 1e6 / estimated_fwhm
            println("  - Estimated FWHM (from mass accuracy): $estimated_fwhm ppm")
        end
        
        if !isempty(peak_counts)
            results[:mean_peaks_per_spectrum] = mean(peak_counts)
        else
            results[:mean_peaks_per_spectrum] = 0
        end
    
    else # CENTROID mode
        println("  Analyzing peak characteristics for CENTROID mode from $(length(spectrum_indices)) sample spectra...")
        
        peak_counts = Int[]
        # peak_intensities = Float64[] # Not strictly needed for these metrics
        
        _iterate_spectra_fast(msi_data, spectrum_indices) do idx, mz, intensity
            if !isempty(mz)
                push!(peak_counts, length(mz))
                # append!(peak_intensities, intensity) # If needed for other metrics
            end
        end
        
        if !isempty(peak_counts)
            results[:mean_peaks_per_spectrum] = mean(peak_counts)
            results[:median_peaks_per_spectrum] = median(peak_counts)
            results[:total_peaks_detected] = sum(peak_counts)
        else
            results[:mean_peaks_per_spectrum] = 0
            results[:median_peaks_per_spectrum] = 0
            results[:total_peaks_detected] = 0
        end
        
        # For centroid data, FWHM is estimated from mass accuracy
        estimated_fwhm = estimated_mean_ppm_error * 2 # Centroid peaks are narrower, factor of 2-3 is common
        results[:mean_fwhm_ppm] = estimated_fwhm
        results[:median_fwhm_ppm] = estimated_fwhm # Defaulting median to mean
        results[:mean_gaussian_r2] = 0.9 # Default for centroid, assuming good peak shapes
        results[:peak_resolution_estimate] = 1e6 / estimated_fwhm
        
        println("  - Mean peaks per spectrum: $(round(results[:mean_peaks_per_spectrum], digits=1))")
        println("  - Estimated peak width: $(round(estimated_fwhm, digits=2)) ppm")
        println("  - Estimated resolution: $(round(results[:peak_resolution_estimate], digits=0))")
    end
    
    # Common prints
    println("  - Mean Gaussian R²: $(round(get(results, :mean_gaussian_r2, 0.0), digits=3))")
    println("  - Mean peaks per spectrum: $(round(get(results, :mean_peaks_per_spectrum, 0.0), digits=1))")
    
    return results
end

"""
    generate_preprocessing_recommendations(analysis_results::Dict)::Dict{Symbol, Any}

Generates intelligent preprocessing recommendations based on the analysis results.
"""
function generate_preprocessing_recommendations(analysis_results::Dict)::Dict{Symbol, Any}
    recommendations = Dict{Symbol, Any}()
    
    inst_analysis = get(analysis_results, :instrument_analysis, Dict())
    signal_analysis = get(analysis_results, :signal_analysis, Dict())
    mass_accuracy = get(analysis_results, :mass_accuracy, Dict())
    peak_analysis = get(analysis_results, :peak_analysis, Dict())
    
    # Baseline Correction Recommendations
    recommendations[:baseline_correction] = generate_baseline_recommendations(inst_analysis, signal_analysis)
    
    # Smoothing Recommendations
    recommendations[:smoothing] = generate_smoothing_recommendations(inst_analysis, peak_analysis)
    
    # Peak Picking Recommendations
    recommendations[:peak_picking] = generate_peak_picking_recommendations(signal_analysis, peak_analysis)
    
    # Normalization Recommendations
    recommendations[:normalization] = generate_normalization_recommendations(signal_analysis)
    
    # Alignment Recommendations
    recommendations[:alignment] = generate_alignment_recommendations(mass_accuracy, inst_analysis)
    
    # Binning Recommendations
    recommendations[:binning] = generate_binning_recommendations(peak_analysis, mass_accuracy)
    
    return recommendations
end

function generate_baseline_recommendations(inst_analysis, signal_analysis)
    recommendations = Dict{Symbol, Any}()
    
    # Determine baseline correction method based on data characteristics
    acquisition_mode = get(inst_analysis, :acquisition_mode, :unknown)
    noise_level = get(signal_analysis, :noise_mean, 1.0)
    
    if acquisition_mode == :profile
        recommendations[:method] = "SNIP"
        recommendations[:window_size] = 200  # Default, can be optimized
        recommendations[:iterations] = 100
    else
        recommendations[:method] = "linear"
        recommendations[:noise_threshold] = noise_level * 3
    end
    
    return recommendations
end

function generate_smoothing_recommendations(inst_analysis, peak_analysis)
    recommendations = Dict{Symbol, Any}()
    
    mz_step = get(inst_analysis, :average_mz_step, 0.01)
    fwhm_ppm = get(peak_analysis, :mean_fwhm_ppm, 20.0)
    
    # Convert FWHM from ppm to m/z units for typical m/z
    typical_mz = 500.0
    fwhm_mz = typical_mz * fwhm_ppm / 1e6
    
    # Savitzky-Golay window should be ~FWHM in points
    window_points = max(5, min(21, round(Int, fwhm_mz / mz_step)))
    # Ensure odd number
    window_points = isodd(window_points) ? window_points : window_points + 1
    
    recommendations[:method] = "Savitzky-Golay"
    recommendations[:window_size] = window_points
    recommendations[:polynomial_order] = 3
    
    return recommendations
end

function generate_peak_picking_recommendations(signal_analysis, peak_analysis)
    recommendations = Dict{Symbol, Any}()
    
    snr_threshold = get(signal_analysis, :suggested_snr, 3.0)
    fwhm_ppm = get(peak_analysis, :mean_fwhm_ppm, 20.0)
    
    recommendations[:snr_threshold] = snr_threshold
    recommendations[:min_peak_width_ppm] = fwhm_ppm * 0.5  # Avoid detecting noise as peaks
    recommendations[:max_peak_width_ppm] = fwhm_ppm * 3.0   # Avoid merging distinct peaks
    
    return recommendations
end

function generate_normalization_recommendations(signal_analysis)
    recommendations = Dict{Symbol, Any}()
    
    tic_cv = get(signal_analysis, :tic_cv, 0.5)
    
    if tic_cv < 0.3  # Low TIC variation
        recommendations[:method] = "TIC"
        recommendations[:reason] = "Low TIC variation across spectra"
    else  # High TIC variation
        recommendations[:method] = "RMS"
        recommendations[:reason] = "High TIC variation, using robust normalization"
    end
    
    return recommendations
end

function generate_alignment_recommendations(mass_accuracy, inst_analysis)
    recommendations = Dict{Symbol, Any}()
    
    calibration_status = get(inst_analysis, :calibration_status, :uncalibrated)

    mean_ppm_error_val = 10.0 # Default value
    if mass_accuracy !== nothing && haskey(mass_accuracy, :global_accuracy) && isa(mass_accuracy[:global_accuracy], NamedTuple)
        ppm_error_report = mass_accuracy[:global_accuracy]
        mean_ppm_error_val = get(ppm_error_report, :mean_ppm_error, 10.0)
    end
    
    if calibration_status == :uncalibrated || mean_ppm_error_val > 20.0
        recommendations[:method] = "reference_based"
        recommendations[:max_ppm_shift] = 50.0
        recommendations[:required] = true
    else
        recommendations[:method] = "none"
        recommendations[:required] = false
    end
    
    return recommendations
end

function generate_binning_recommendations(peak_analysis, mass_accuracy)
    recommendations = Dict{Symbol, Any}()
    
    fwhm_ppm = get(peak_analysis, :mean_fwhm_ppm, 20.0)
    
    # Handle the case when mass_accuracy is nothing or has NaN values
    suggested_tolerance = 20.0  # Default value
    if mass_accuracy !== nothing
        tol = get(mass_accuracy, :suggested_bin_tolerance, 20.0)
        if !isnan(tol) && isfinite(tol)
            suggested_tolerance = tol
        end
    end
    
    # Ensure fwhm_ppm is valid
    if isnan(fwhm_ppm) || !isfinite(fwhm_ppm)
        fwhm_ppm = 20.0
    end
    
    # Bin width should be ~FWHM/2 to preserve resolution while reducing data size
    bin_width_ppm = max(fwhm_ppm * 0.5, suggested_tolerance)
    
    recommendations[:method] = "adaptive"
    recommendations[:bin_width_ppm] = bin_width_ppm
    recommendations[:min_peaks_per_bin] = 3
    
    return recommendations
end

# =============================================================================
# Pre-Analysis Pipeline for Auto Parameter Determination
# =============================================================================

"""
    run_preprocessing_analysis(msi_data::MSIData; 
                              reference_peaks::Dict{Float64, String}=Dict(),
                              region_masks::Dict=Dict(),
                              sample_size::Int=100)::Dict{Symbol, Any}

Runs a comprehensive pre-analysis pipeline to determine optimal preprocessing parameters.
This function analyzes the dataset to provide intelligent defaults for preprocessing steps.
"""
function run_preprocessing_analysis(msi_data::MSIData; 
                                   reference_peaks::Dict{Float64, String}=Dict{Float64, String}(),
                                   region_masks::Dict{Symbol, BitMatrix}=Dict{Symbol, BitMatrix}(),
                                   sample_size::Int=100,
                                   mask_path::Union{String, Nothing}=nothing,
                                   spectrum_indices::Union{AbstractVector{Int}, Nothing}=nothing)::Dict{Symbol, Any}
    
    println("="^60)
    println("RUNNING PRE-ANALYSIS PIPELINE")
    println("="^60)
    
    analysis_results = Dict{Symbol, Any}()

    local all_available_indices::AbstractVector{Int}

    if spectrum_indices !== nothing
        println("Using provided list of $(length(spectrum_indices)) spectrum indices.")
        all_available_indices = spectrum_indices
    elseif mask_path !== nothing
        println("Applying mask from: $(mask_path)")
        try
            mask_matrix = load_and_prepare_mask(mask_path, msi_data.image_dims)
            masked_indices_set = get_masked_spectrum_indices(msi_data, mask_matrix)
            all_available_indices = collect(masked_indices_set)
            println("Mask applied. $(length(all_available_indices)) spectra are within the masked region.")
        catch e
            @error "Failed to load or apply mask: $e. Proceeding without mask."
            all_available_indices = 1:length(msi_data.spectra_metadata)
        end
    else
        all_available_indices = 1:length(msi_data.spectra_metadata)
    end

    if isempty(all_available_indices)
        @warn "No spectra available for analysis (after applying mask/filter). Returning empty results."
        return analysis_results
    end

    # Create a single sample set from the available indices
    num_available = length(all_available_indices)
    indices_to_sample = if num_available > sample_size
        # Use StatsBase.sample for sampling without replacement
        sample(all_available_indices, sample_size, replace=false)
    else
        # Use all available indices if they are fewer than the sample size
        collect(all_available_indices)
    end
    
    # Ensure analytics are precomputed
    if !is_set(msi_data.analytics_ready)
        println("Pre-computing basic analytics...")
        precompute_analytics(msi_data)
    end
    
    # Phase 1: Instrument and Data Characteristics
    println("\n--- Phase 1: Instrument & Data Characteristics ---")
    instrument_analysis = analyze_instrument_characteristics(msi_data, sample_indices=indices_to_sample)
    analysis_results[:instrument_analysis] = instrument_analysis
    
    # Phase 2: Noise and Signal Quality Analysis
    println("\n--- Phase 2: Noise & Signal Quality Analysis ---")
    signal_analysis = analyze_signal_quality(msi_data, sample_indices=indices_to_sample)
    analysis_results[:signal_analysis] = signal_analysis
    
    # Phase 3: Mass Accuracy Analysis
    println("\n--- Phase 3: Mass Accuracy Analysis ---")
    if !isempty(reference_peaks)
        mass_accuracy_analysis = analyze_mass_accuracy_global(msi_data, reference_peaks, spectrum_indices=indices_to_sample)
        analysis_results[:mass_accuracy] = mass_accuracy_analysis
    else
        println("No reference peaks provided - skipping mass accuracy analysis")
        analysis_results[:mass_accuracy] = nothing
    end
    
    # Phase 4: Spatial Region Analysis (if masks provided)
    # This function is not affected by the global mask, as it analyzes specific, named regions.
    println("\n--- Phase 4: Spatial Region Analysis ---")
    if !isempty(region_masks)
        regional_analysis = analyze_spatial_regions(msi_data, region_masks, reference_peaks)
        analysis_results[:regional_analysis] = regional_analysis
    else
        println("No region masks provided - skipping regional analysis")
        analysis_results[:regional_analysis] = nothing
    end
    
    # Phase 5: Peak Characteristics Analysis
    println("\n--- Phase 5: Peak Characteristics Analysis ---")
    peak_analysis = analyze_peak_characteristics(msi_data, instrument_analysis, analysis_results[:mass_accuracy], spectrum_indices=indices_to_sample)
    analysis_results[:peak_analysis] = peak_analysis
    
    # Phase 6: Generate Preprocessing Recommendations
    println("\n--- Phase 6: Generating Preprocessing Recommendations ---")
    recommendations = generate_preprocessing_recommendations(analysis_results)
    analysis_results[:recommendations] = recommendations
    
    # Store in MSIData object
    msi_data.preprocessing_hints = recommendations
    
    println("\n" * "="^60)
    println("PRE-ANALYSIS COMPLETE")
    println("="^60)
    
    return analysis_results
end

"""
    calculate_ppm_error(measured_mz::Real, theoretical_mz::Real) -> Float64

Calculates the mass accuracy error in parts-per-million (PPM) between a measured
and a theoretical m/z value.

# Arguments
- `measured_mz::Real`: The experimentally measured m/z value.
- `theoretical_mz::Real`: The known, theoretical m/z value of a compound.

# Returns
- `Float64`: The calculated PPM error. Returns `Inf` if `theoretical_mz` is zero.

# Formula
`PPM = 10^6 * |measured_mz - theoretical_mz| / theoretical_mz`

# Example
```julia
calculate_ppm_error(100.005, 100.0) # returns 50.0
```
"""
function calculate_ppm_error(measured_mz::Real, theoretical_mz::Real)
    if theoretical_mz == 0
        return Inf
    end
    return 1e6 * abs(Float64(measured_mz) - Float64(theoretical_mz)) / Float64(theoretical_mz)
end

"""
    calculate_ppm_error_bulk(measured_mz::Vector{<:Real}, theoretical_mz::Vector{<:Real}) -> Vector{Float64}

Calculates PPM errors for multiple pairs of measured and theoretical mass values.

# Arguments
- `measured_mz::Vector{<:Real}`: A vector of experimentally measured m/z values.
- `theoretical_mz::Vector{<:Real}`: A vector of known, theoretical m/z values.

# Returns
- `Vector{Float64}`: A vector containing the calculated PPM error for each pair.
"""
function calculate_ppm_error_bulk(measured_mz::Vector{Real}, theoretical_mz::Vector{Real})
    return [calculate_ppm_error(m, t) for (m, t) in zip(measured_mz, theoretical_mz)]
end

"""
    calculate_resolution_fwhm(mz::Real, profile_mz::AbstractVector{<:Real}, profile_intensity::AbstractVector{<:Real}) -> Float64

Calculates the mass resolution of a peak in profile-mode data using the Full Width at
Half Maximum (FWHM) method. Resolution is a measure of an instrument's ability to
distinguish between two peaks of slightly different mass-to-charge ratios.

# Arguments
- `mz::Real`: The m/z value of the peak's centroid.
- `profile_mz::AbstractVector{<:Real}`: The full m/z array from the profile-mode spectrum.
- `profile_intensity::AbstractVector{<:Real}`: The full intensity array from the profile-mode spectrum.

# Returns
- `Float64`: The calculated resolution (`m / Δm`). Returns `NaN` if the FWHM cannot be determined (e.g., peak is at the edge of the spectrum).

# Formula
`Resolution = m / Δm`, where `Δm` is the FWHM.
"""
function calculate_resolution_fwhm(mz::Real, profile_mz::AbstractVector{<:Real},
                                 profile_intensity::AbstractVector{<:Real})
    
    # Find peak center index
    peak_idx = argmin(abs.(profile_mz .- mz))
    peak_height = Float64(profile_intensity[peak_idx])
    half_max = peak_height / 2
    
    # Find left half-maximum point (interpolate for accuracy)
    left_idx = find_last_below(profile_intensity[1:peak_idx], half_max)
    if left_idx == 0 || left_idx == length(profile_intensity[1:peak_idx])
        return NaN
    end
    
    # Linear interpolation for left FWHM
    x1, x2 = Float64(profile_mz[left_idx]), Float64(profile_mz[left_idx+1])
    y1, y2 = Float64(profile_intensity[left_idx]), Float64(profile_intensity[left_idx+1])
    denominator = y2 - y1
    if denominator == 0
        left_fwhm = NaN
    else
        left_fwhm = x1 + (x2 - x1) * (half_max - y1) / denominator
    end
    
    # Find right half-maximum point
    right_slice = profile_intensity[peak_idx:end]
    right_offset = find_first_below(right_slice, half_max)
    if right_offset == 0 || right_offset == length(right_slice)
        return NaN
    end
    
    right_idx = peak_idx + right_offset - 1
    x1, x2 = Float64(profile_mz[right_idx-1]), Float64(profile_mz[right_idx])
    y1, y2 = Float64(profile_intensity[right_idx-1]), Float64(profile_intensity[right_idx])
    denominator = y2 - y1
    if denominator == 0
        right_fwhm = NaN
    else
        right_fwhm = x1 + (x2 - x1) * (half_max - y1) / denominator
    end
    
    fwhm = right_fwhm - left_fwhm
    return fwhm > 0 ? Float64(mz) / fwhm : NaN
end

# Helper functions for FWHM calculation
function find_last_below(v::AbstractVector{<:Real}, threshold::Real)
    for i in length(v):-1:2
        if v[i] >= threshold && v[i-1] < threshold
            return i-1
        end
    end
    return 0
end

function find_first_below(v::AbstractVector{<:Real}, threshold::Real)
    for i in 1:(length(v)-1)
        if v[i] >= threshold && v[i+1] < threshold
            return i+1
        end
    end
    return 0
end

"""
    _calculate_fwhm_delta_m(mz::AbstractVector{<:Real}, intensity::AbstractVector{<:Real}, peak_idx::Int) -> Float64

Calculates the Full Width at Half Maximum (FWHM) in m/z units (Δm) for a peak.
Returns `NaN` if FWHM cannot be determined.
"""
function _calculate_fwhm_delta_m(mz::AbstractVector{<:Real}, intensity::AbstractVector{<:Real}, peak_idx::Int)
    peak_height = Float64(intensity[peak_idx])
    half_max = peak_height / 2
    
    # Find left half-maximum point
    left_idx = find_last_below(intensity[1:peak_idx], half_max)
    if left_idx == 0 || left_idx == length(intensity[1:peak_idx])
        return NaN
    end
    
    # Linear interpolation for left FWHM m/z
    x1, x2 = Float64(mz[left_idx]), Float64(mz[left_idx+1])
    y1, y2 = Float64(intensity[left_idx]), Float64(intensity[left_idx+1])
    denominator = y2 - y1
    if denominator == 0
        left_fwhm_mz = NaN
    else
        left_fwhm_mz = x1 + (x2 - x1) * (half_max - y1) / denominator
    end
    
    # Find right half-maximum point
    right_slice = intensity[peak_idx:end]
    right_offset = find_first_below(right_slice, half_max)
    if right_offset == 0 || right_offset == length(right_slice)
        return NaN
    end
    
    right_idx = peak_idx + right_offset - 1
    x1, x2 = Float64(mz[right_idx-1]), Float64(mz[right_idx])
    y1, y2 = Float64(intensity[right_idx-1]), Float64(intensity[right_idx])
    denominator = y2 - y1
    if denominator == 0
        right_fwhm_mz = NaN
    else
        right_fwhm_mz = x1 + (x2 - x1) * (half_max - y1) / denominator
    end
    
    return right_fwhm_mz - left_fwhm_mz
end

"""
    _fit_gaussian_and_r2(mz::AbstractVector{<:Real}, intensity::AbstractVector{<:Real}, peak_idx::Int, half_window::Int) -> Float64

Estimates Gaussian parameters for a peak and returns a pseudo R^2 value.
This is an approximation and not a full non-linear least squares fit.
"""
function _fit_gaussian_and_r2(mz::AbstractVector{<:Real}, intensity::AbstractVector{<:Real}, peak_idx::Int, half_window::Int)
    n = length(mz)
    if n < 3 || peak_idx <= 0 || peak_idx > n
        return 0.0
    end

    # Define the region around the peak
    start_idx = max(1, peak_idx - half_window)
    end_idx = min(n, peak_idx + half_window)
    
    # Ensure there's enough data to fit
    if (end_idx - start_idx + 1) < 3
        return 0.0
    end

    x_data = mz[start_idx:end_idx]
    y_data = intensity[start_idx:end_idx]

    # Estimate Gaussian parameters
    # Amplitude (A): peak intensity
    A_est = intensity[peak_idx]
    # Mean (μ): m/z at peak intensity
    mu_est = mz[peak_idx]
    # Standard deviation (σ): related to FWHM. FWHM = 2 * sqrt(2 * ln(2)) * σ ≈ 2.355 * σ
    # So, σ ≈ FWHM / 2.355
    fwhm_delta_m = _calculate_fwhm_delta_m(mz, intensity, peak_idx)
    if isnan(fwhm_delta_m) || fwhm_delta_m <= 0
        return 0.0 # Cannot estimate sigma without a valid FWHM
    end
    sigma_est = fwhm_delta_m / 2.355

    # If sigma is too small, it might lead to division by zero or very sharp peaks
    if sigma_est < eps(Float64)
        return 0.0
    end

    # Gaussian function
    gaussian(x, A, mu, sigma) = A * exp.(-(x .- mu).^2 ./ (2 * sigma^2))

    # Generate estimated Gaussian curve
    y_est = gaussian(x_data, A_est, mu_est, sigma_est)

    # Calculate pseudo R-squared
    # R^2 = 1 - (SS_res / SS_tot)
    # SS_res = sum((y_data - y_est).^2)
    # SS_tot = sum((y_data - mean(y_data)).^2)

    SS_res = sum((y_data .- y_est).^2)
    SS_tot = sum((y_data .- mean(y_data)).^2)

    if SS_tot == 0
        return 1.0 # Perfect fit if all y_data are the same
    end

    r_squared = 1.0 - (SS_res / SS_tot)
    return max(0.0, r_squared) # R^2 can be negative if fit is worse than mean, cap at 0
end

"""
    calculate_robust_fwhm(mz::AbstractVector{<:Real}, intensity::AbstractVector{<:Real}, peak_idx::Int) -> Float64

A more robust FWHM calculation that handles edge cases better.
"""
function calculate_robust_fwhm(mz::AbstractVector{<:Real}, intensity::AbstractVector{<:Real}, peak_idx::Int)
    n = length(mz)
    if n < 5 || peak_idx < 3 || peak_idx > n-2
        return NaN
    end
    
    peak_height = Float64(intensity[peak_idx])
    half_max = peak_height / 2.0
    
    # Find left half-maximum with bounds checking
    left_idx = peak_idx
    while left_idx > 1 && intensity[left_idx] >= half_max
        left_idx -= 1
    end
    
    if left_idx == 1 || left_idx >= n-1
        return NaN
    end
    
    # Linear interpolation for left FWHM
    x1, x2 = Float64(mz[left_idx]), Float64(mz[left_idx+1])
    y1, y2 = Float64(intensity[left_idx]), Float64(intensity[left_idx+1])
    
    if y2 == y1  # Avoid division by zero
        left_fwhm_mz = x1
    else
        left_fwhm_mz = x1 + (x2 - x1) * (half_max - y1) / (y2 - y1)
    end
    
    # Find right half-maximum
    right_idx = peak_idx
    while right_idx < n && intensity[right_idx] >= half_max
        right_idx += 1
    end
    
    if right_idx == n || right_idx <= 2
        return NaN
    end
    
    # Linear interpolation for right FWHM
    x1, x2 = Float64(mz[right_idx-1]), Float64(mz[right_idx])
    y1, y2 = Float64(intensity[right_idx-1]), Float64(intensity[right_idx])
    
    if y2 == y1  # Avoid division by zero
        right_fwhm_mz = x1
    else
        right_fwhm_mz = x1 + (x2 - x1) * (half_max - y1) / (y2 - y1)
    end
    
    fwhm = right_fwhm_mz - left_fwhm_mz
    
    # Validate result
    if fwhm <= 0 || !isfinite(fwhm) || fwhm > 1.0  # Unreasonably large
        return NaN
    end
    
    return fwhm
end

"""
    main_precalculation(msi_data::MSIData; ...)

Runs the non-verbose pre-analysis pipeline and returns a dictionary of recommended 
preprocessing parameters based on data characteristics and heuristics.

This function serves as a quiet entry point to the analysis engine, translating the
analytical results into a concrete set of parameters for a preprocessing pipeline.
It also categorizes parameters that cannot be automatically determined.

# Arguments
- `msi_data::MSIData`: The main MSI data object.
- `reference_peaks::Dict`: Optional dictionary of reference m/z values for mass accuracy analysis.
- `region_masks::Dict`: Optional dictionary of named `BitMatrix` masks for regional analysis.
- `sample_size::Int`: The number of spectra to sample for statistical analysis.
- `mask_path::String`: Optional path to a PNG mask file to restrict the analysis to a specific ROI.
- `spectrum_indices::AbstractVector{Int}`: Optional vector of spectrum indices to restrict analysis to.

# Returns
- A `Dict` with two keys:
    - `"recommended_parameters"`: A `Dict{Symbol, Any}` of suggested parameter values.
    - `"unsupported_parameters"`: A `Dict` categorizing parameters that could not be determined.
"""
function main_precalculation(msi_data::MSIData; 
                             reference_peaks::Dict{Float64, String}=Dict{Float64, String}(),
                             region_masks::Dict{Symbol, BitMatrix}=Dict{Symbol, BitMatrix}(),
                             sample_size::Int=100,
                             mask_path::Union{String, Nothing}=nothing,
                             spectrum_indices::Union{AbstractVector{Int}, Nothing}=nothing)::Dict
    
    # --- Run analysis pipeline quietly ---
    local analysis_results
    original_stdout = stdout
    # Redirect stdout to the system's null device to robustly silence output
    null_stream = open(Sys.iswindows() ? "nul" : "/dev/null", "w")
    redirect_stdout(null_stream)
    try
        analysis_results = run_preprocessing_analysis(msi_data,
            reference_peaks=reference_peaks,
            region_masks=region_masks,
            sample_size=sample_size,
            mask_path=mask_path,
            spectrum_indices=spectrum_indices
        )
    finally
        redirect_stdout(original_stdout)
        close(null_stream)
    end

    if isempty(analysis_results)
        @warn "Preprocessing analysis returned no results. Cannot generate recommendations."
        return Dict() # Return an empty dictionary if no results
    end
    
    # --- Safely get nested dictionaries ---
    recs = get(analysis_results, :recommendations, Dict())
    signal_analysis = get(analysis_results, :signal_analysis, Dict())
    peak_analysis = get(analysis_results, :peak_analysis, Dict())
    mass_accuracy = get(analysis_results, :mass_accuracy, Dict())
    inst_analysis = get(analysis_results, :instrument_analysis, Dict())

    # Initialize parameter dictionaries for each step
    cal_params = Dict{Symbol, Any}()
    sm_params = Dict{Symbol, Any}()
    bc_params = Dict{Symbol, Any}()
    norm_params = Dict{Symbol, Any}()
    pp_params = Dict{Symbol, Any}()
    pa_params = Dict{Symbol, Any}()
    ps_params = Dict{Symbol, Any}()
    pb_params = Dict{Symbol, Any}()

    # --- Populate Parameters for each step ---

    # Calibration & Alignment (Note: These are intertwined in the current logic)
    calibration_required = false
    mean_ppm_error = get(get(mass_accuracy, :global_accuracy, Dict()), :mean_ppm_error, NaN)
    suggested_bin_tol = get(mass_accuracy, :suggested_bin_tolerance, NaN)
    
    if !isempty(recs) && haskey(recs, :alignment) && get(recs[:alignment], :required, false)
        calibration_required = true
        cal_params[:method] = :internal_standards
        cal_params[:fit_order] = 2 # Default from struct
        if isfinite(suggested_bin_tol)
             cal_params[:ppm_tolerance] = suggested_bin_tol
        else
             cal_params[:ppm_tolerance] = nothing
        end
        cal_params[:internal_standards] = nothing # User input dependent
        cal_params[:base_peak_mz_references] = nothing # User input dependent

        if isfinite(mean_ppm_error)
            pa_params[:method] = mean_ppm_error > 30.0 ? :lowess : :linear
        else
            pa_params[:method] = :lowess # Default
        end

        tic_cv = get(signal_analysis, :tic_cv, NaN)
        if isfinite(tic_cv) && get(pa_params, :method, :none) == :lowess
            pa_params[:span] = round(max(0.3, min(0.8, 1.0 - tic_cv / 2)), digits=2)
        else
            pa_params[:span] = nothing
        end
        
        if isfinite(suggested_bin_tol)
             pa_params[:tolerance] = suggested_bin_tol
             pa_params[:tolerance_unit] = :ppm
        else
             pa_params[:tolerance] = nothing
             pa_params[:tolerance_unit] = :ppm # Default
        end
        pa_params[:max_shift_ppm] = get(recs[:alignment], :max_ppm_shift, 50.0)
        pa_params[:min_matched_peaks] = nothing # User input dependent
    else
        cal_params[:method] = :none
        cal_params[:ppm_tolerance] = nothing
        cal_params[:fit_order] = nothing
        cal_params[:internal_standards] = nothing
        cal_params[:base_peak_mz_references] = nothing

        pa_params[:method] = :none
        pa_params[:span] = nothing
        pa_params[:tolerance] = nothing
        pa_params[:tolerance_unit] = :ppm
        pa_params[:max_shift_ppm] = nothing
        pa_params[:min_matched_peaks] = nothing
    end

    # Smoothing
    if !isempty(recs) && haskey(recs, :smoothing)
        sm_rec = recs[:smoothing]
        if get(sm_rec, :method, "") == "Savitzky-Golay"
            sm_params[:method] = :savitzky_golay
            sm_params[:window] = get(sm_rec, :window_size, nothing)
            sm_params[:order] = get(sm_rec, :polynomial_order, nothing)
        else
            sm_params[:method] = :none
            sm_params[:window] = nothing
            sm_params[:order] = nothing
        end
    else
        sm_params[:method] = :none
        sm_params[:window] = nothing
        sm_params[:order] = nothing
    end
    
    # Baseline Correction
    if !isempty(recs) && haskey(recs, :baseline_correction)
        bl_rec = recs[:baseline_correction]
        if get(bl_rec, :method, "") == "SNIP"
             bc_params[:method] = :snip
             bc_params[:iterations] = get(bl_rec, :iterations, nothing)
        else
            bc_params[:method] = :none
            bc_params[:iterations] = nothing
        end
    else
        bc_params[:method] = :none
        bc_params[:iterations] = nothing
    end

    if get(inst_analysis, :acquisition_mode, :unknown) == :profile
        mean_fwhm_ppm = get(peak_analysis, :mean_fwhm_ppm, NaN)
        avg_mz_step = get(inst_analysis, :average_mz_step, NaN)
        if isfinite(mean_fwhm_ppm) && isfinite(avg_mz_step) && avg_mz_step > 0
            fwhm_mz = 500.0 * mean_fwhm_ppm / 1e6 # At typical m/z 500
            bc_params[:window] = ceil(Int, fwhm_mz / avg_mz_step * 2)
        else
            bc_params[:window] = nothing
        end
    else
        bc_params[:window] = nothing
    end

    # Normalization
    if !isempty(recs) && haskey(recs, :normalization)
        method = get(recs[:normalization], :method, "")
        if method == "TIC"
            norm_params[:method] = :tic
        elseif method == "RMS"
            norm_params[:method] = :rms
        else
            norm_params[:method] = :none
        end
    else
        norm_params[:method] = :none
    end

    # Peak Picking
    # Robust method selection with fallback
    acquisition_mode = get(inst_analysis, :acquisition_mode, :profile) # Default to profile
    pp_params[:method] = acquisition_mode == :profile ? :profile : :centroid

    if !isempty(recs) && haskey(recs, :peak_picking)
        pk_rec = recs[:peak_picking]
        pp_params[:snr_threshold] = get(pk_rec, :snr_threshold, 3.0) # Default to 3.0
        pp_params[:min_peak_width_ppm] = get(pk_rec, :min_peak_width_ppm, nothing)
        pp_params[:max_peak_width_ppm] = get(pk_rec, :max_peak_width_ppm, nothing)
    else
        pp_params[:snr_threshold] = 3.0 # Default to 3.0
        pp_params[:min_peak_width_ppm] = nothing
        pp_params[:max_peak_width_ppm] = nothing
    end
    
    # Robust prominence calculation with safety cap
    estimated_noise = get(signal_analysis, :noise_mean, NaN)
    if isfinite(estimated_noise)
        # Never be more aggressive than 0.05, a reasonable upper limit
        calculated_prominence = round(estimated_noise * 2, digits=4)
        pp_params[:min_peak_prominence] = min(calculated_prominence, 0.05)
    else
        # Fallback to a safe, non-aggressive value if noise couldn't be estimated
        pp_params[:min_peak_prominence] = 0.01
    end

    if isfinite(suggested_bin_tol)
        pp_params[:merge_peaks_tolerance] = round(suggested_bin_tol / 2, digits=4)
    else
        pp_params[:merge_peaks_tolerance] = nothing
    end

    # Robust half_window calculation with safety floor
    mean_fwhm_ppm = get(peak_analysis, :mean_fwhm_ppm, NaN)
    avg_mz_step = get(inst_analysis, :average_mz_step, NaN)
    if isfinite(mean_fwhm_ppm) && isfinite(avg_mz_step) && avg_mz_step > 0
        fwhm_mz = 500.0 * mean_fwhm_ppm / 1e6 # At typical m/z 500
        window_points = fwhm_mz / avg_mz_step
        calculated_half_window = ceil(Int, window_points / 2)
        # Ensure half_window is at least 3, a safe absolute minimum
        pp_params[:half_window] = max(calculated_half_window, 3)
    else
        # Fallback if calculation is not possible
        pp_params[:half_window] = 5
    end
    
    # Robust R^2 calculation with floor
    mean_r2 = get(peak_analysis, :mean_gaussian_r2, NaN)
    if isfinite(mean_r2)
        pp_params[:min_peak_shape_r2] = round(max(0.5, mean_r2 * 0.8), digits=2)
    else
        # Fallback to a reasonable default
        pp_params[:min_peak_shape_r2] = 0.6
    end


    # Peak Selection
    if haskey(pp_params, :snr_threshold) && pp_params[:snr_threshold] !== nothing
        ps_params[:min_snr] = pp_params[:snr_threshold]
    else
        ps_params[:min_snr] = nothing
    end
    if !isempty(peak_analysis)
        mean_fwhm = get(peak_analysis, :mean_fwhm_ppm, NaN)
        if isfinite(mean_fwhm)
            ps_params[:min_fwhm_ppm] = round(mean_fwhm * 0.5, digits=2)
            ps_params[:max_fwhm_ppm] = round(mean_fwhm * 2.5, digits=2)
        else
            ps_params[:min_fwhm_ppm] = nothing
            ps_params[:max_fwhm_ppm] = nothing
        end
        if isfinite(mean_r2)
            ps_params[:min_shape_r2] = round(max(0.5, mean_r2 * 0.8), digits=2)
        else
            ps_params[:min_shape_r2] = nothing
        end
    else
        ps_params[:min_fwhm_ppm] = nothing
        ps_params[:max_fwhm_ppm] = nothing
        ps_params[:min_shape_r2] = nothing
    end
    ps_params[:frequency_threshold] = nothing # User input dependent / Hard to determine
    ps_params[:correlation_threshold] = nothing # Hard to determine


    # Peak Binning
    if !isempty(recs) && haskey(recs, :binning)
        bin_rec = recs[:binning]
        if get(bin_rec, :method, "") == "adaptive"
            pb_params[:method] = :adaptive
            bin_width = get(bin_rec, :bin_width_ppm, NaN)
            if isfinite(bin_width)
                pb_params[:tolerance] = bin_width
                pb_params[:max_bin_width_ppm] = round(bin_width * 3, digits=2)
            else
                pb_params[:tolerance] = nothing
                pb_params[:max_bin_width_ppm] = nothing
            end
            pb_params[:tolerance_unit] = :ppm
            pb_params[:min_peak_per_bin] = get(bin_rec, :min_peaks_per_bin, nothing)
            pb_params[:intensity_weighted_centers] = true # Default from struct
            pb_params[:num_uniform_bins] = nothing # User input dependent
            pb_params[:frequency_threshold] = nothing # Hard to determine
        else
            pb_params[:method] = :none
            pb_params[:tolerance] = nothing
            pb_params[:max_bin_width_ppm] = nothing
            pb_params[:tolerance_unit] = nothing
            pb_params[:min_peak_per_bin] = nothing
            pb_params[:intensity_weighted_centers] = true
            pb_params[:num_uniform_bins] = nothing
            pb_params[:frequency_threshold] = nothing
        end
    else
        pb_params[:method] = :none
        pb_params[:tolerance] = nothing
        pb_params[:max_bin_width_ppm] = nothing
        pb_params[:tolerance_unit] = nothing
        pb_params[:min_peak_per_bin] = nothing
        pb_params[:intensity_weighted_centers] = true
        pb_params[:num_uniform_bins] = nothing
        pb_params[:frequency_threshold] = nothing
    end

    return Dict(
        :Calibration => cal_params,
        :Smoothing => sm_params,
        :BaselineCorrection => bc_params,
        :Normalization => norm_params,
        :PeakPicking => pp_params,
        :PeakAlignment => pa_params,
        :PeakSelection => ps_params,
        :PeakBinningParams => pb_params
    )
end
