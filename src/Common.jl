# src/Common.jl

# --- Unified Error Types ---
abstract type MSIError <: Exception end

struct FileFormatError <: MSIError
    msg::String
end

struct SpectrumNotFoundError <: MSIError
    id
end

struct InvalidSpectrumError <: MSIError
    id
    reason::String
end

# --- Shared Constants ---
const DEFAULT_CACHE_SIZE = 100
const DEFAULT_NUM_BINS = 2000


"""
    validate_spectrum_data(mz, intensity, id)

Checks a spectrum for common data integrity issues.

- Checks for `NaN` or `Inf` values in intensity and m/z arrays.
- Verifies that the m/z array is sorted in ascending order.
- Ensures that m/z and intensity arrays have the same length.

# Throws
- `InvalidSpectrumError` if any of the checks fail.
"""
function validate_spectrum_data(mz::AbstractVector, intensity::AbstractVector, id)
    if length(mz) != length(intensity)
        throw(InvalidSpectrumError(id, "m/z and intensity arrays have different lengths"))
    end

    if !all(isfinite, mz) || !all(isfinite, intensity)
        throw(InvalidSpectrumError(id, "Spectrum contains NaN or Inf values"))
    end

    if !issorted(mz)
        throw(InvalidSpectrumError(id, "m/z array is not sorted"))
    end

    return true
end
