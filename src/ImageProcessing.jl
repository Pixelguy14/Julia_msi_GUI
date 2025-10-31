module ImageProcessing

using Images
using ImageBinarization
using ImageMorphology
using ImageComponentAnalysis

export process_image_pipeline

# ===================================================================
# CORE PROCESSING PIPELINE
# ===================================================================

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

end
