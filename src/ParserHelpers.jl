# src/ParserHelpers.jl

"""
This file provides common helper functions for parsing mzML and imzML files.
"""

# ============================================================================
#
# Data Structures and Low-Level Parser Helpers
#
# ============================================================================

"""
    SpecDim

A struct to hold configuration for a spectral data axis (e.g., m/z or intensity).
"""
mutable struct SpecDim
    Format::Type
    Packed::Bool
    Axis::Int
    Skip::Int
end

"""
    find_tag(stream, regex::Regex)

Reads a stream line-by-line until a line matches the provided regex.

# Returns
- A `RegexMatch` object if a match is found, otherwise throws an error.
"""
function find_tag(stream, regex::Regex)
    while !eof(stream)
        line = readline(stream)
        isTag = match(regex, line)
        if isTag !== nothing
            return isTag
        end
    end
    error("Tag not found for regex: $regex")
end

"""
    get_attribute(source::String, tag::String = "([^=]+)")

Retrieves an attribute's value from an XML tag string.

# Returns
- A `RegexMatch` object containing the attribute and its value.
"""
function get_attribute(source::AbstractString, tag::String = "([^=]+)")
    # Construct the regex pattern string
    pattern_str = "\\s" * tag * "=\"([^\"]*)\""
    regStr = Regex(pattern_str)
    return match(regStr, source)
end

"""
    configure_spec_dim(stream)

Parses a block of `<cvParam>` tags from an XML stream to configure a `SpecDim`
struct. It reads accessions to determine the data format (e.g., `Float32`),
compression status (`zlib`), and axis type (m/z vs. intensity).

# Arguments
- `stream`: An IO stream positioned at the start of the `cvParam` block.

# Returns
- A `SpecDim` struct populated with the parsed configuration.
"""
function configure_spec_dim(stream)
    axis = SpecDim(Float64, false, 1, 0)
    offset = position(stream)

    while !eof(stream)
        currLine = readline(stream)
        matchInfo = match(r"^\s*<(cvParam)", currLine)

        if matchInfo === nothing
            matchInfo = match(r"^\s*", currLine)
            axis.Skip = position(stream) - offset - length(currLine) + length(matchInfo.match)
            return axis
        end

        index = length(matchInfo.captures[1])
        attr_match = get_attribute(currLine[index:end], "accession")

        if attr_match !== nothing
            accession = attr_match.captures[1]
            if accession == "MS:1000515"       # intensity array
                axis.Axis = 2
            elseif accession == "MS:1000519"   # 32-bit integer
                axis.Format = Int32
            elseif accession == "MS:1000521"   # 32-bit float
                axis.Format = Float32
            elseif accession == "MS:1000522"   # 64-bit integer
                axis.Format = Int64
            elseif accession == "MS:1000574"   # zlib compression
                axis.Packed = true
            end
        end
    end
    return axis # Should be unreachable if file is well-formed
end

# ============================================================================
#
# Data Structures and Helpers for mzML lazy parsing
#
# ============================================================================

"""
    CVParams

A mutable struct used during the parsing of mzML files to temporarily hold
CV parameter information for a data array before it is converted into a
`SpectrumAsset`.

# Fields
- `format`: The data type (e.g., `Float64`).
- `is_compressed`: A boolean indicating if the data is compressed.
- `axis_type`: A symbol indicating the array type (`:mz` or `:intensity`).
"""
mutable struct CVParams
    format::Type
    is_compressed::Bool
    axis_type::Symbol
end
