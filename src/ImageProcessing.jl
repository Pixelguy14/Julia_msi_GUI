# src/ImageProcessing

using Images
using ImageBinarization
using ImageMorphology
using ImageComponentAnalysis
using Colors # For converting to grayscale

export process_image_pipeline
export load_and_prepare_mask # Export the new function

"""
    load_and_prepare_mask(mask_path::String, target_dims::Tuple{Int, Int})

Loads a PNG image mask, converts it to a binary (Boolean) matrix,
and resizes it to the specified `target_dims`. White pixels in the mask
are considered `true` (part of the ROI), and black pixels are `false`.

# Arguments
- `mask_path`: Absolute path to the PNG mask file.
- `target_dims`: A tuple `(width, height)` representing the desired output dimensions.

# Returns
- A `BitMatrix` of `target_dims` where `true` indicates the ROI.
"""
function load_and_prepare_mask(mask_path::String, target_dims::Tuple{Int, Int})
    if !isfile(mask_path)
        error("Mask file not found: $(mask_path)")
    end

    # Load the image
    img = Images.load(mask_path)

    # Convert to grayscale if it's a color image
    gray_img = Gray.(img)

    # Binarize using Otsu's method - white regions become 'true'
    binary_img = binarize(gray_img, Otsu())

    # Resize to target dimensions
    resized_img = imresize(binary_img, (target_dims[2], target_dims[1]))

    # Ensure the output is a BitMatrix, as imresize can change the type
    return resized_img .> 0.5
end

# ===================================================================
# CORE PROCESSING PIPELINE
# ===================================================================

"""
    process_image_pipeline(gray_img; otsu_scale=1.0, noise_size_percent=0.1, hole_size_percent=0.05, smoothing=2)

Applies a multi-step image processing pipeline to a grayscale image to segment regions of interest.

The pipeline consists of:
1.  **Binarization**: An adjusted Otsu's threshold is used to create a binary image.
2.  **Noise Removal**: Small white regions (noise) are removed using an area opening operation.
3.  **Hole Filling**: Small black regions (holes) within larger objects are filled.
4.  **Edge Smoothing**: The edges of the final regions are smoothed using a morphological closing operation.

# Arguments
- `gray_img`: The input grayscale image (`Matrix{<:Gray}`).

# Keyword Arguments
- `otsu_scale`: A factor to scale the automatically determined Otsu threshold. Values > 1.0 make the threshold stricter (less white), < 1.0 make it more lenient (more white). Default: `1.0`.
- `noise_size_percent`: The percentage of the total image area used as a threshold to remove small white noise components. Default: `0.1`.
- `hole_size_percent`: The percentage of the total image area used as a threshold to fill black holes in white components. Default: `0.05`.
- `smoothing`: The size of the kernel for the final edge smoothing (closing) operation. Default: `2`.

# Returns
- A tuple containing four images representing the intermediate steps of the pipeline:
    1.  `binary_img`: The result of the initial binarization.
    2.  `noise_removed_img`: The image after noise removal.
    3.  `holes_filled_img`: The image after filling holes.
    4.  `smoothed_img`: The final smoothed image.
"""
function process_image_pipeline(gray_img;
                                otsu_scale=1.0,
                                noise_size_percent=0.1,
                                hole_size_percent=0.05,
                                smoothing=2)
    
    # --- Step 1: Otsu Binarization ---
    otsu_threshold = find_threshold(gray_img, Otsu())
    adjusted_threshold = otsu_threshold * otsu_scale
    binary_img = gray_img .>= adjusted_threshold
    
    # --- Smart Parameter Scaling ---
    image_area = length(gray_img)
    noise_size_pixels = round(Int, image_area * noise_size_percent)
    hole_size_pixels = round(Int, image_area * hole_size_percent)

    # --- Step 2: Remove Small White Regions (Noise) ---
    # area_opening is the correct morphological operation for this.
    noise_removed_img = area_opening(binary_img, min_area=noise_size_pixels)

    # --- Step 3: Fill Small Black Holes ---
    # area_closing is the dual of area_opening and fills holes.
    holes_filled_img = area_closing(noise_removed_img, min_area=hole_size_pixels)

    # --- Step 4: Smooth Edges ---
    # A morphological closing with a small disk smooths outlines.
    smoothing_kernel = ones(Bool, (smoothing, smoothing))
    smoothed_img = closing(holes_filled_img, smoothing_kernel)

    # --- Return all intermediate steps for visualization ---
    return binary_img, noise_removed_img, holes_filled_img, smoothed_img
end
