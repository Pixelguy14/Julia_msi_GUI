# A Comprehensive Guide to Mass Spectrometry Imaging Preprocessing
## Rationale, Methods, and Implementation

---

### Table of Contents
1. [Introduction](#introduction)
2. [Profile vs. Centroid Mode Data](#profile-vs-centroid-mode-data)
3. [General Principles of MSI Preprocessing](#general-principles-of-msi-preprocessing)
4. [Step-by-Step Breakdown of Preprocessing Steps](#step-by-step-breakdown-of-preprocessing-steps)
    * [Validation and Quality Control](#validation-and-quality-control)
    * [Intensity Transformation (Stabilization)](#intensity-transformation-stabilization)
    * [Smoothing](#smoothing)
    * [Baseline Correction](#baseline-correction)
    * [Normalization](#normalization)
    * [Peak Detection](#peak-detection)
    * [Peak Selection (Filtering)](#peak-selection-filtering)
    * [Calibration](#calibration)
    * [Peak Alignment](#peak-alignment)
    * [Peak Binning and Feature Matrix Generation](#peak-binning-and-feature-matrix-generation)
5. [Example Walkthrough: From Raw Spectra to Feature Matrix](#example-walkthrough-from-raw-spectra-to-feature-matrix)
    * [Phase 1 – Instrument & Data Characteristics](#phase-1--instrument--data-characteristics)
    * [Phase 2 – Noise & Signal Quality](#phase-2--noise--signal-quality)
    * [Phase 3 – Mass Accuracy Analysis](#phase-3--mass-accuracy-analysis)
    * [Phase 4 – Spatial Region Analysis](#phase-4--spatial-region-analysis)
    * [Phase 5 – Peak Characteristics](#phase-5--peak-characteristics)
    * [Phase 6 – Preprocessing Recommendations](#phase-6--preprocessing-recommendations)
6. [Extra Tips for the Analyst](#extra-tips-for-the-analyst)
    * [Isotopic and Charge State Matters](#isotopic-and-charge-state-matters)
    * [Always Visualize the Average Spectrum First](#always-visualize-the-average-spectrum-first)
    * [Preprocessing Is Not a One-Shot Process](#preprocessing-is-not-a-one-shot-process)

---

## Introduction
Mass Spectrometry Imaging (MSI) generates thousands of mass spectra, each acquired from a distinct pixel location. In the literature, before any interpretation can be drawn from these data, a series of preprocessing steps must be applied. These steps correct for instrumental artifacts, reduce noise, and align data across pixels to produce a clean, comparable feature matrix suitable for statistical analysis and image reconstruction.

Inside the **JuliaMSI** code (`Preprocessing.jl` and `PreprocessingPipeline.jl`), we have a preprocessing workflow inspired by established R packages like **MALDIquant** and **Cardinal**.

This document explains the rationale behind each step, the required inputs, the scientific meaning of parameters, and the context behind key algorithms such as **TIC normalization**, **PQN**, **wavelet-based peak detection**, and **LOWESS alignment**.

---

## Profile vs. Centroid Mode Data

Mass spectrometry data can be acquired in two primary modes: **profile** and **centroid**. Understanding the difference is crucial because preprocessing steps must be adapted accordingly.

*   **Profile mode**:
    The instrument records intensity at every digitized point along the *m/z* axis, producing a continuous waveform. Peaks appear as Gaussian-like shapes. This mode retains full peak shape information, which is valuable for accurate mass measurement, resolution assessment, and advanced peak detection. However, file sizes are large and noise is present between peaks.
    
*   **Centroid mode**:
    The instrument performs real-time peak detection and reduces each peak to a single centroid (*m/z*, intensity) pair, discarding the underlying profile. This dramatically reduces file size but loses peak shape details. Centroid data are essentially a list of detected peaks with no information about peak width or baseline.

In **profile mode**, algorithms must operate on the full *m/z* and intensity vectors, performing tasks like smoothing, baseline correction, and peak picking on continuous data. In **centroid mode**, these steps are often unnecessary or simplified: smoothing and baseline correction are not applicable because the data are already peak lists; peak picking is replaced by filtering based on intensity thresholds. 

> [!IMPORTANT]
> The JuliaMSI code internally distinguishes this difference in some steps by acquiring the mode from the metadata on the files, but it's up to the analyst to decide skipping smoothing/baseline steps for centroid data, for example.

---

## General Principles of MSI Preprocessing

All preprocessing steps share a common goal: **to end up with an aligned feature matrix where rows correspond to pixels (spectra) and columns correspond to *m/z* bins (features).**

The steps are applied sequentially, often in a specific order to avoid compounding errors. For example, **baseline correction should precede normalization** because baseline offsets can distort total ion current calculations. Similarly, **peak alignment should follow peak detection** because alignment requires a set of common peaks across spectra.

---

## Step-by-Step Breakdown of Preprocessing Steps

### Validation and Quality Control
*(This step doesn't require analyst input)*
*   **Functions in code**: `validate_spectrum`, `qc_is_empty`, `qc_is_regular`
*   **Purpose**: Ensure that each input spectrum is well-formed and contains usable data before any processing.
*   **Inputs required**:
    *   `mz::AbstractVector{<:Real}`: The *m/z* values.
    *   `intensity::AbstractVector{<:Real}`: The corresponding intensities.

> **Why these inputs?**
> Every subsequent algorithm depends on the integrity of these vectors. Checking for non-emptiness, equal length, finite values, non-negativity, and monotonic *m/z* order prevents cryptic errors later.

---

### Intensity Transformation (Stabilization)
*   **Function**: `transform_intensity_core`
*   **Purpose**: Apply a variance-stabilizing transformation to make the data more homoscedastic (constant variance across intensity range).
*   **Inputs required**:
    *   `intensity::AbstractVector{<:Real}`
    *   `method::Symbol` (e.g., `:sqrt`, `:log1p`)

> **Why these inputs?**
> The transformation operates on the intensity vector only; *m/z* values remain unchanged. The choice of method depends on the data characteristics (e.g., presence of zeros).
>
> Stabilization uses functions like **square root**, **log**, **log2**, **log10** and **log1p** helps the analyst stabilize variance for count data, or handle quadratical variance. The analyst must know which one of these methods to use according to the data they are visualizing (e.g., by the mean spectrum plot, or a plot from a specific coordinate).
>
> In MSI in general, **square root stabilization** is common because ion counting statistics approximate Poisson.
>
> **Example**: A spectrum with intense peaks and low baseline noise: after square root, the dynamic range is compressed, making low-intensity features more comparable to high-intensity ones.

---

### Smoothing
*   **Function**: `smooth_spectrum_core`, `moving_average_smooth`
*   **Purpose**: Reduce high-frequency noise while preserving peak shapes.
*   **Inputs required**:
    *   `y::AbstractVector{<:Real}`: Intensity vector.
    *   `method::Symbol` (`:savitzky_golay`, `:moving_average`)
    *   `window::Int`: Size of the smoothing window.
    *   `order::Int`: Polynomial order for Savitzky–Golay.

> **Why these inputs?**
> Smoothing operates on the intensity vector only. The window size must be an odd integer for symmetry; the polynomial order must be less than the window. These constraints ensure the filter behaves properly (e.g., Savitzky–Golay requires the window bigger than the order parameter).
>
> **Savitzky–Golay** fits a low-degree polynomial to the local window, preserving moments of the peak (height, width) better than a simple moving average, which can distort peak shapes.
>
> **Example**: A noisy peak with 5% random noise: after Savitzky–Golay, the peak becomes smoother while its centroid and height remain accurate.

---

### Baseline Correction
*   **Function**: `apply_baseline_correction_core`, `_snip_baseline_impl`, `convex_hull_baseline`, `median_baseline`
*   **Purpose**: Estimate and subtract the background signal (baseline) from each spectrum.
*   **Inputs required**:
    *   `y::AbstractVector{<:Real}`: Intensity vector.
    *   `method::Symbol` (`:snip`, `:convex_hull`, `:median`)
    *   `iterations::Int` (for SNIP)
    *   `window::Int` (for median)

> **Why these inputs?**
> Baseline correction requires only the intensity vector. The method-specific parameters control the aggressiveness and shape of the baseline estimate.
>
> *   **SNIP (Sensitive Nonlinear Iterative Peak clipping)**: Iteratively replaces each point with the minimum of itself and the average of its neighbors. After many iterations, only the baseline remains. SNIP is the default in **MALDIquant** due to its effectiveness.
> *   **Convex hull**: Constructs the upper convex hull of the spectrum (or lower hull for baseline) and interpolates.
> *   **Moving median**: Estimates baseline as the local median. Simple but can be biased by dense peaks.
>
> **Example**: In a spectrum with a rising baseline after SNIP correction, the baseline is removed, revealing peaks sitting on a flat background.

---

### Normalization
*   **Function**: `apply_normalization_core`, `tic_normalize`, `pqn_normalize`, `median_normalize`, `rms_normalize`
*   **Purpose**: Correct for differences in total ion abundance between spectra, making intensities comparable across pixels.
*   **Inputs required**:
    *   `y::AbstractVector{<:Real}` (per-spectrum) or a matrix `M` (for PQN).
    *   `method::Symbol`

> **Why these inputs?**
> Normalization operates on intensities. **PQN** requires a matrix of spectra (rows = *m/z* bins, columns = spectra) because it computes a reference spectrum across all samples.
>
> *   **TIC (Total Ion Current)**: Divides each intensity by the sum of all intensities in the spectrum. A simple algorithm but sensitive to a few very intense peaks (e.g., matrix peaks).
> *   **Median normalization**: Divides by the median intensity. Robust to outliers but can be unstable if many intensities are zero.
> *   **RMS normalization**: Divides by root-mean-square, preserving the shape of the spectrum while scaling to unit energy.
> *   **PQN (Probabilistic Quotient Normalization)**: Calculates a reference spectrum and scales each spectrum by the median of the quotients (intensity/reference). This method is preferred when some peaks vary drastically because it downweights their influence.
>
> **Example**: Two spectra from different tissue regions: one has a huge lipid peak, the other has lower overall signal. TIC would suppress the lipid peak in the first spectrum, possibly losing biological information. PQN would scale both spectra more equitably.

---

### Peak Detection
*   **Functions**: `detect_peaks_profile_core`, `detect_peaks_wavelet_core`, `detect_peaks_centroid_core`
*   **Purpose**: Identify peaks (signals of interest) in each spectrum.
*   **Inputs required**:
    *   *m/z*, intensity vectors.
    *   Method-specific parameters: `snr_threshold`, `half_window`, `min_peak_prominence`, `merge_peaks_tolerance`, etc.

> **Why these inputs?**
> Peak detection uses both *m/z* and intensity to locate local maxima and evaluate peak quality. Parameters control sensitivity and filtering. A peak is a local maximum that rises significantly above the local noise. Detection algorithms vary by data mode:
>
> *   **Profile mode (local maxima)**: Scans for points that are highest in a sliding window (`half_window`). Then filters by **SNR** (signal-to-noise ratio) and prominence (height relative to surrounding minima). Close peaks are merged if within `merge_peaks_tolerance`.
> *   **Wavelet-based**: Uses **Continuous Wavelet Transform (CWT)** to identify peaks at multiple scales, making it robust to varying peak widths. Local maxima in the CWT magnitude indicate peak positions. Wavelet is often used for complex spectra with overlapping peaks. CWT-based peak detection is implemented in **MALDIquant** as the default method for profile data.
> *   **Centroid mode**: Simply filters the existing peak list by SNR, as peaks are already identified.
>
> **Example**: A spectrum with a small peak on the shoulder of a large peak: wavelet detection may resolve both, while simple local maxima might miss the shoulder.

---

### Peak Selection (Filtering)
*   **Function**: `apply_peak_selection`
*   **Purpose**: Remove peaks that do not meet quality criteria, such as minimum SNR, acceptable width, or good Gaussian shape.
*   **Inputs required**:
    *   peaks list (from peak detection).
    *   Thresholds: `min_snr`, `min_fwhm_ppm`, `max_fwhm_ppm`, `min_shape_r2`.

> **Why these inputs?**
> These thresholds reflect the analyst's confidence in peak quality. They are typically derived from instrument specifications or empirical data.
>
> Not all local maxima are true peaks. Noise spikes, baseline artifacts, or electronic interference can produce spurious peaks. Filtering by **SNR** removes noise. **FWHM** (full width at half maximum) in **ppm** ensures that peaks have physically reasonable widths (e.g., isotopic peaks have predictable spacing). **Shape R²** from a Gaussian fit assesses how well the peak resembles a theoretical peak; poor shape may indicate saturation or interference.
>
> **Example**: A peak with SNR < 3 is likely noise and is discarded. A peak with FWHM > 100 *ppm* might be background rather than a specific ion.

---

### Calibration
*(This step requires that you have a list of Internal Standards inside JuliaMSI preprocessing tab filled in)*
*   **Function**: `calibrate_spectra_core`, `find_calibration_peaks_core`
*   **Purpose**: Correct systematic mass errors by aligning the *m/z* axis to known reference masses (internal standards).
*   **Inputs required**:
    *   spectra (list of `MutableSpectrum`).
    *   `internal_standards`: dictionary of theoretical *m/z* values.
    *   `ppm_tolerance`: matching tolerance.
    *   `fit_order`: polynomial order for calibration curve.

> **Why these inputs?**
> Calibration requires a set of known masses that are present in the spectra. The algorithm matches detected peaks to these references, then fits a transformation (e.g., linear) to map measured to theoretical *m/z*.
>
> Mass spectrometers drift; this step is made to correct this drift, ensuring that the same *m/z* value corresponds to the same ion across the entire imaging dataset. **Linear or quadratic calibration** is standard. The code currently uses linear interpolation between matched peaks, which is sufficient for most MSI applications where drift is small.
>
> **Example**: A reference peak expected at **760.585 Da** (a lipid) is found at **760.590 Da** in a spectrum. Calibration shifts the entire axis so that the peak aligns with the theoretical value.

---

### Peak Alignment
*   **Function**: `align_peaks_lowess_core`
*   **Purpose**: Correct residual *m/z* shifts between spectra after calibration, ensuring that the same analyte peak appears at the same *m/z* across all pixels.
*   **Inputs required**:
    *   Reference spectrum's peak *m/z* list.
    *   Target spectrum's peak *m/z* list.
    *   Parameters: method (`:lowess`, `:linear`, `:ransac`), `tolerance`, `tolerance_unit`, `max_shift_ppm`.

> **Why these inputs?**
> Alignment operates on peak lists, not full spectra, to be efficient. It finds common peaks between a reference spectrum (e.g., the one with most peaks) and each target, then fits a warping function.
>
> Even after calibration, small nonlinear shifts may persist due to sample topography or matrix effects. Alignment brings all spectra onto a common *m/z* scale, which is essential for constructing a feature matrix.
>
> **LOWESS (Locally Weighted Scatterplot Smoothing)**: Fits a smooth curve through the matched peak pairs, allowing for nonlinear corrections. It is robust to outliers and does not assume a global polynomial shape. It has been used in **MALDIquant**. **RANSAC (Random Sample Consensus)** is a robust alternative. The code defaults to LOWESS because it is well-suited for smooth, continuous drift.
>
> **Example**: Peaks in a target spectrum are systematically shifted by **+0.01 Da** at low *m/z* and **+0.005 Da** at high *m/z*. LOWESS fits a curve that corrects each peak individually, producing aligned *m/z* values.

---

### Peak Binning and Feature Matrix Generation
*   **Function**: `bin_peaks_core`
*   **Purpose**: Group peaks from all spectra into common *m/z* bins, producing a final feature matrix.
*   **Inputs required**:
    *   spectra (with detected peaks).
    *   `params::PeakBinning`: method (`:adaptive` or `:uniform`), `tolerance`, `frequency_threshold`, etc.

> **Why these inputs?**
> Binning requires the list of peaks from each spectrum. The method determines how bins are defined.
>
> *   **Adaptive binning**: Clusters peaks based on their *m/z* proximity (using tolerance). Bins are created where peaks are dense, and bins without enough peaks are discarded.
> *   **Uniform binning**: Divides the total *m/z* range into a fixed number of equally spaced bins. 
>
> After alignment, peaks from different spectra still have slightly different *m/z* values due to residual noise. Binning assigns them to discrete bins, creating a matrix suitable for multivariate analysis. The **frequency threshold** removes bins that appear in too few spectra (likely noise).
>
> **Example**: Peaks around *m/z* **760.58** from 100 spectra are grouped into one bin. The intensity for each spectrum in that bin is the maximum intensity of peaks assigned to it.

---

## Example Walkthrough: From Raw Spectra to Feature Matrix
### with Intelligent Parameter Recommendations

To illustrate the entire pipeline, consider a mass spectrometry imaging dataset of a mouse brain section acquired in **profile mode**. The analyst has a list of internal standards for lipids (e.g., known phospholipids at *m/z* **760.585**, **798.540**, **834.528**) to aid in mass calibration.

Before running the full preprocessing, the **JuliaMSI** code executes the intelligent pre-analysis pipeline (`main_precalculation` which is a wrapper for `run_preprocessing_analysis`). This function samples a representative subset of spectra (e.g., 100 pixels) and performs a battery of diagnostic tests to suggest optimal parameters.

#### The pre-analysis yields a report:

*   **Phase 1 – Instrument & Data Characteristics**:
    *   Acquisition mode: **profile** (continuous waveforms)
    *   *m/z* axis type: **regular** (constant step of ~0.01 Da)
    *   Dynamic range: ~4 orders of magnitude
*   **Phase 2 – Noise & Signal Quality**:
    *   Estimated noise level (median absolute deviation): 45 counts
    *   Suggested SNR threshold: **3.0**
    *   TIC coefficient of variation: **35%** → indicates moderate heterogeneity, so robust normalization (e.g., **PQN**) may be beneficial.
*   **Phase 3 – Mass Accuracy Analysis** (using the provided lipid standards):
    *   Mean PPM error: **12.3 ppm**
    *   Standard deviation: **8.1 ppm**
    *   Suggested bin tolerance (mean + 3σ): **36.6 ppm** (capped between 10 and 100 ppm).
*   **Phase 4 – Spatial Region Analysis**:
    *   *(if region masks were supplied) – skipped in this example.*
*   **Phase 5 – Peak Characteristics**:
    *   Mean FWHM: **22.5 ppm** (converted from Δm at *m/z* 500)
    *   Mean Gaussian R²: **0.85** (peaks are well-shaped)
    *   Mean peaks per spectrum: **187**
*   **Phase 6 – Preprocessing Recommendations**:
    *   The system synthesizes all information into a dictionary of suggested parameters for each step. Parameters that cannot be determined automatically (e.g., frequency threshold for binning, minimum matched peaks) are left as `nothing`; the analyst must supply them based on experimental goals.

> [!NOTE]
> If you have internal standards, remember to re-run the preprocessing recommendations after you have loaded them. This is only for precalculation and is not a required step for analysts that know about their data.

#### Execution of the Pipeline:
The analyst can now execute the pipeline, possibly overriding a few parameters (e.g., increasing SNR threshold to 4 for more stringent peak picking).

1.  **Validation**: All spectra pass QC; one empty spectrum is flagged and removed.
2.  **Intensity transformation (sqrt)**: Stabilizes variance, making low-intensity features more comparable.
3.  **Baseline correction (SNIP, 100 iterations)**: Removes background using the recommended window of 200 points.
4.  **Smoothing (Savitzky–Golay, window 9, order 2)**: Reduces noise while preserving peak shape.
5.  **Peak detection (profile mode)**: Detects peaks with SNR ≥ 3, prominence ≥ 0.002, and uses the recommended half-window of 5 points. Peaks closer than 18 *ppm* are merged.
6.  **Peak selection**: Filters peaks: keeps only those with SNR ≥ 3, FWHM between 11 and 55 *ppm*, and shape R² ≥ 0.68.
7.  **Calibration**: Using the internal standards, each spectrum’s *m/z* axis is linearly corrected based on matched peaks within 36.6 *ppm* tolerance.
8.  **Peak alignment (LOWESS)**: Aligns all spectra to a reference spectrum (the one with most peaks) with a span of 0.6, using the same 36.6 *ppm* tolerance. This corrects residual nonlinear drift.
9.  **Normalization (RMS)**: Divides each spectrum by its root-mean-square intensity, chosen because TIC variation was high.
10. **Peak binning (adaptive)**: Groups peaks within 36.6 *ppm* into bins. Bins containing fewer than 3 peaks or appearing in <20% of spectra are discarded. The final bin centers are intensity-weighted averages.
11. **Feature matrix**: A **999 × 425 matrix** is produced, where rows are pixels and columns are *m/z* bins. This matrix is ready for PCA, segmentation, or statistical testing.

#### How the Analyst Uses the Recommendations:
*   The suggested values are data-driven, but the analyst can override them based on domain knowledge. Always verify; experimenting is encouraged as data integrity is kept.
*   Parameters like `frequency_threshold` in peak selection and binning require the analyst’s judgment. A threshold of **0.2** means that a bin must contain a peak in at least 20% of all pixels to be retained. This reduces noise but may discard rare but real signals.
*   This is an iterative refinement; the analyst may rerun the analysis with adjusted parameters (e.g., a lower bin tolerance).
*   JuliaMSI aims for a tailored preprocessing workflow that respects both the instrument’s characteristics and the biological question.

---

## Extra Tips for the Analyst

The preprocessing pipeline is powerful, but its success ultimately depends on the analyst’s informed choices. Here are practical tips to guide you through the process, especially when working with internal standards and interpreting intermediate results.

> [!WARNING]
> This is a really resource-intensive pipeline dependent on the sample, so run your precautions.

### Isotopic and Charge State Matters
When providing a dictionary of reference masses for calibration, ensure that the theoretical *m/z* values correspond to the same charge state and isotopic composition as the peaks you expect to see in your data.
*   **Example**: If your instrument detects lipids mainly as **[M+H]⁺** ions, the reference mass must be the monoisotopic mass plus a proton. Using the neutral mass will introduce a systematic offset of **1.0078 Da**, leading to poor calibration. (Make the same adjustment for negative ion mode).
*   If isotopic peaks are prominent, you may want to include multiple isotopes (e.g., **M**, **M+1**, **M+2**) as separate references, but be aware that their relative intensities vary. The calibration algorithm will match the closest detected peak. This is common for MSI preprocessing.
*   The pre-analysis suggests a bin tolerance based on observed mass error. A tolerance that is too large may match the wrong peak (e.g., an isobaric interference); a tolerance that is too small may miss valid matches.

### Always Visualize the Average Spectrum First
Before running any preprocessing, generate a sum or average spectrum. This single plot tells you:
1.  **Baseline shape**: If it's flat, rising, or undulating, it can guide the baseline correction method (**SNIP** works well for most shapes; **convex hull** is faster but may fail on curvy baselines).
2.  **Noise level**: By zooming in.
3.  ***m/z* Scatter**: How scattered are the spectra *m/z*.
4.  **Peak density & Intensity range**: Influence the peak detection thresholds (SNR, prominence) and binning strategy.

*   If the average spectrum shows a **rising baseline** from low to high *m/z*, SNIP with ~100 iterations will handle it well. If extremely wavy, increase iterations or try a different method.
*   If running from the console (`run_preprocessing.jl`), the output of `run_preprocessing_analysis` provides:
    *   **TIC coefficient of variation (CV)**: A CV > 30% suggests substantial pixel-to-pixel variation, favoring **RMS** or **PQM** over TIC.
    *   **Mean PPM error**: If >20 *ppm*, calibration is essential. If <5 *ppm* and drift is linear, you might skip alignment.
    *   **FWHM and Gaussian R²**: Indicate peak quality. If R² < 0.6, peaks may be poorly shaped; consider increasing the `min_peak_shape_r2` threshold.
*   **Acquisition mode**: For **centroid data**, never apply smoothing or baseline correction because those steps assume continuous waveforms.
    *   The pipeline will set `:method = :centroid` in peak picking and use permissive filters.
    *   Use a low SNR threshold (e.g., 2.0).
    *   Disable shape-based filtering (`min_peak_shape_r2 = 0.0`).
    *   Centroid data do not retain peak width information; set `min_fwhm_ppm = 0.0` and `max_fwhm_ppm = 500.0` to avoid discarding real peaks.
    *   Adaptive binning still works; tolerance should be based on expected mass accuracy.

### Preprocessing Is Not a One-Shot Process
The recommended parameters are a starting point, not a final answer.
*   **Reproducibility**: Keep a record of all parameters used. The `params` dictionary can be saved as a **JSON**. This is essential for publications.
*   **Biological Context**: If looking for rare metabolites appearing in few pixels, lower the `frequency_threshold` in binning. If interested only in high-abundance lipids, increase the SNR threshold.
*   **Modularity**: The pipeline is modular; you can always rerun parts of it with adjusted parameters without starting from scratch.

---
*Any revision to the code, bug, error, report to julian.sierrag@icloud.com, if you we're to cite this file, remember i provide no citation and its based on internet knowledge of the web about these algorithms, so do your own research and cite properly. to cite JuliaMSI platform, please cite:*

José Julián Sierra-Álvarez, Martín Orlando Camargo-Escalante, Carlos Daniel Sierra-Álvarez, Carmelo Hernández-Caricio, Juan Francisco Moreno-Luna, Isabel Buendía-Corona, Robert Winkler,
JuliaMSI: A high-performance graphical platform for mass spectrometry imaging data analysis,
Analytica Chimica Acta,
Volume 1377,
2025,
344613,
ISSN 0003-2670,
https://doi.org/10.1016/j.aca.2025.344613
