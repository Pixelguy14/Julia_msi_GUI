# src/ParserHelpers.jl

"""
This file provides common helper functions for parsing mzML and imzML files.
"""

# ============================================================================
#
# Data Structures and Low-Level Parser Helpers
#
# ============================================================================

# Pre-compiled Regex patterns for efficiency
const ACCESSION_REGEX = r"accession=\"([^\"]*)\""
const VALUE_REGEX = r"value=\"([^\"]*)\""
const NAME_REGEX = r"name=\"([^\"]*)\""
const CVP_START_REGEX = r"^\s*<cvParam"
const CVP_END_TAG_REGEX = r"^\s*</"
const SPECTRUM_TAG_ID_REGEX = r"<spectrum\s+[^>]*id=\"([^\"]+)\"" # For direct ID extraction
const BINARY_DATA_ARRAY_ENCODED_LENGTH_REGEX = r"<binaryDataArray\s+[^>]*encodedLength=\"(\d+)\""
const INDEX_SPECTRUM_TAG_REGEX = r"<index\s+name=\"spectrum\""
const OFFSET_TAG_REGEX = r"<offset[^>]*>(\d+)</offset>"
const INDEX_LIST_OFFSET_REGEX = r"<indexListOffset>(\d+)</indexListOffset>"
const COUNT_REGEX = r"count=\"(\d+)\""
const REF_PARAM_GROUP_LIST_REGEX = r"<referenceableParamGroupList"
const REF_PARAM_GROUP_ID_REGEX = r"<referenceableParamGroup id=\"([^\"]+)\""
const X_COORD_REGEX = r"IMS:1000050[^>]*value=\"(\d+)\""
const Y_COORD_REGEX = r"IMS:1000051[^>]*value=\"(\d+)\""
const DEFAULT_ARRAY_LENGTH_REGEX = r"defaultArrayLength=\"(\d+)\""
const ENCODED_LENGTH_ATTR_REGEX = r"encodedLength=\"(\d+)\""
const ARRAY_LENGTH_CV_REGEX = r"IMS:1000103.*?value=\"(\d+)\""
const ENCODED_LEN_CV_REGEX = r"IMS:1000104.*?value=\"(\d+)\""
const OFFSET_CV_REGEX = r"IMS:1000102.*?value=\"(\d+)\""


"""
    SpecDim

A struct to hold configuration for a spectral data axis (e.g., m/z or intensity).
"""
mutable struct SpecDim
    Format::Type
    Packed::Bool
    Axis::Int
    Skip::Int
    Mode::SpectrumMode
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
    axis = SpecDim(Float64, false, 1, 0, UNKNOWN)
    
    while !eof(stream)
        currLine = readline(stream)
        
        # Check for start of cvParam
        matchInfo = match(CVP_START_REGEX, currLine)

        if matchInfo === nothing
            # Check for end of cvParam block (any closing tag)
            if match(CVP_END_TAG_REGEX, currLine) !== nothing
                return axis
            end
            continue
        end

        # Extract accession using the pre-compiled regex
        attr_match = match(ACCESSION_REGEX, currLine)

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
            elseif accession == "MS:1000127"   # centroid spectrum
                axis.Mode = CENTROID
            elseif accession == "MS:1000128"   # profile spectrum  
                axis.Mode = PROFILE
            end
        end
    end
    return axis
end
