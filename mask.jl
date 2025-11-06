module MaskApp

# --- Dependencies ---
using GenieFramework
using Images, ImageBinarization, ImageMorphology, ImageComponentAnalysis
using NativeFileDialog, FileIO, ImageCore, Printf, Dates, JSON
using Pkg, Libz, PlotlyBase, CairoMakie, Colors, Base64
using Statistics, NaturalSort, LinearAlgebra, StipplePlotly
using Base.Filesystem: mv

using MSI_src
using .MSI_src: MSIData, process_image_pipeline

# Plot Handling
include("./julia_imzML_visual.jl")

# Image Processing Pipeline
using ImageBinarization

function load_and_binarize_mask(path)
    img = load(path)
    gray_img = ensure_grayscale(img)
    # Binarize to ensure only pure black and white values
    return binarize(gray_img, Otsu())
end

function ensure_grayscale(img)
    if eltype(img) <: AbstractRGB; return Gray.(img); end
    if eltype(img) <: AbstractRGBA; return Gray.(RGB.(img)); end
    if eltype(img) <: Color; return Gray.(img); end
    return img
end

function alter_image(img_path, otsu_scale, noise_size_percent, hole_size_percent, smoothing_level)
    try
        # Strip query parameters from img_path
        clean_img_path = replace(img_path, r"\?.*" => "")
        full_image_path = joinpath("public", lstrip(clean_img_path, '/'))
        if !isfile(full_image_path)
            return (success=false, message="Image file not found: $(full_image_path)", path="")
        end

        original_img = load(full_image_path)
        gray_img = Float32.(ensure_grayscale(original_img))

        binary, noise_removed, holes_filled, smoothed =
            process_image_pipeline(gray_img;
                otsu_scale=otsu_scale, noise_size_percent=noise_size_percent,
                hole_size_percent=hole_size_percent, smoothing=smoothing_level)

        output_dir = joinpath("public", "css", "masks")
        mkpath(output_dir)

        path_smooth = "/css/masks/smoothed.png"
        save(joinpath(output_dir, basename(path_smooth)), smoothed)
        
        return (success=true, message="Image processed successfully!", path=path_smooth)

    catch e
        @error "Mask editor processing failed" exception=(e, catch_backtrace())
        return (success=false, message="Error processing image: $(sprint(showerror, e))", path="")
    end
end

# New helper function to display a single slice
function display_slice(slice_path)
    if isempty(slice_path)
        return [PlotlyBase.heatmap(x=Vector{Float64}(), y=Vector{Float64}())], PlotlyBase.Layout(margin=attr(l=0,r=0,t=0,b=0,pad=0)), false
    end
    
    full_path = joinpath("public", lstrip(slice_path, '/'))
    if !isfile(full_path)
        return [PlotlyBase.heatmap(x=Vector{Float64}(), y=Vector{Float64}())], PlotlyBase.Layout(margin=attr(l=0,r=0,t=0,b=0,pad=0)), false
    end

    plotdata, plotlayout, _, _ = loadImgPlot(slice_path, "", 0.0) # No mask overlay
    return plotdata, plotlayout, true
end

function get_timestamped_path(base_path)
    clean_path = replace(base_path, r"\?.*" => "")
    return clean_path * "?t=" * string(time_ns())
end

function refresh_editor_preview(imgInt, smoothed_mask_path, imgTrans)
    if !isempty(smoothed_mask_path)
        timestamped_path = get_timestamped_path(smoothed_mask_path)
        plotdata, plotlayout, _, _ = loadImgPlot(imgInt, timestamped_path, imgTrans)
        return plotdata, plotlayout
    else
        return [PlotlyBase.heatmap(x=Vector{Float64}(), y=Vector{Float64}())], PlotlyBase.Layout(margin=attr(l=0,r=0,t=0,b=0,pad=0))
    end
end

function flood_fill!(img::AbstractMatrix, x::Int, y::Int, fill_color)
    h, w = size(img)
    if !(1 <= x <= h && 1 <= y <= w); return; end
    
    target_color = img[x, y]
    if target_color == fill_color; return; end
    
    q = [(x, y)]
    img[x, y] = fill_color
    
    while !isempty(q)
        cx, cy = popfirst!(q)
        
        for (dx, dy) in [(0, 1), (0, -1), (1, 0), (-1, 0)]
            nx, ny = cx + dx, cy + dy
            if 1 <= nx <= h && 1 <= ny <= w && img[nx, ny] == target_color
                img[nx, ny] = fill_color
                push!(q, (nx, ny))
            end
        end
    end
end

function update_main_plot(imgInt::String, smoothed_mask_path::String, imgTrans::Float64)
    try
        if !isempty(imgInt) && !isempty(smoothed_mask_path)
            # Both slice and mask available - show overlay
            plotdata, plotlayout, _, _ = loadImgPlot(imgInt, smoothed_mask_path, imgTrans)
            return plotdata, plotlayout, true
        elseif !isempty(imgInt)
            # Only slice available - show slice alone
            plotdata, plotlayout, show = display_slice(imgInt)
            return plotdata, plotlayout, show
        elseif !isempty(smoothed_mask_path)
            # Only mask available - show mask alone
            plotdata, plotlayout, _, _ = loadImgPlot(smoothed_mask_path, "", 0.0)
            return plotdata, plotlayout, true
        else
            # Nothing to show
            return [PlotlyBase.heatmap(x=Vector{Float64}(), y=Vector{Float64}())], 
                   PlotlyBase.Layout(margin=attr(l=0,r=0,t=0,b=0,pad=0)), 
                   false
        end
    catch e
        @error "Failed to update main plot" exception=(e, catch_backtrace())
        return [PlotlyBase.heatmap(x=Vector{Float64}(), y=Vector{Float64}())], 
               PlotlyBase.Layout(margin=attr(l=0,r=0,t=0,b=0,pad=0)), 
               false
    end
end

@genietools

@app begin
    # --- State for Slice Selection (from app.jl) ---
    @in selected_folder_main = ""
    @out available_folders = String[]
    @out image_available_folders = String[]
    @private registry_init_done = false

    @out imgInt = "" # Path to current slice
    @out current_msi = ""
    @out msgimg = "Please select a dataset."
    @in imgTrans = 0.5 # Default transparency
    
    # --- State for Mask Workflow ---
    @in otsu_scale = 1.0
    @in noise_size_percent = 0.01
    @in hole_size_percent = 0.005
    @in smoothing_level = 3
    @in is_editing_mask = false
    @in is_browsing_slices = true

    @in mask_input_path = "" # Reactive trigger for processing
    @out smoothed_mask_path = "" # Path to the final smoothed.png

    @out plotdata_verify = [PlotlyBase.heatmap(x=Vector{Float64}(), y=Vector{Float64}())]
    @out plotlayout_verify = PlotlyBase.Layout(margin=attr(l=0,r=0,t=0,b=0,pad=0))
    @out show_verification_plot = false

    @out show_editor = false # To launch the custom editor

    # --- Messages and Warnings ---
    @out mask_editor_message = ""
    @out mask_editor_warning = false
    @out progress = false

    # --- Buttons ---
    @in btn_use_slice_as_mask = false
    @in btn_upload_mask = false
    @in btn_edit_manually = false
    @in btn_img_plus = false
    @in btn_img_minus = false
    @in btn_change_slice = false # New button
    @in btn_flip_mask = false
    @in btn_save_final_mask = false
    @in btn_move_mask = false

    # --- Manual Editor State ---
    @in brush_size = 10
    @in brush_color = "#ff0000"  # Red brush
    @in current_tool = "brush"   # "brush", "eraser", "bucket", "drag"
    @in editor_scale = 1.0
    @in bucket_fill_trigger = Dict()
    
    # Canvas dimensions matching the smoothed mask
    @out canvas_width = 512
    @out canvas_height = 512
    
    # Mouse coordinates for drawing
    @in mouse_x = 0
    @in mouse_y = 0
    @in is_drawing = false
    
    # Editor operations
    @in rotate_degrees = 0
    @in flip_direction = "horizontal"
    @in move_direction = "right"
    @in move_pixels = 10
    @in move_mask_payload = Dict()
    
    # Final save
    @in final_mask_name = ""
    @in updated_mask_data = "" # New property for canvas data
    @in canvas_refresh_trigger = 0 # New property to trigger canvas refresh

    # --- Handlers ---

    @onchange updated_mask_data begin
        if !isempty(updated_mask_data)
            try
                # Data URL format: data:image/png;base64,iVBORw0KGgo...
                # Extract base64 part
                base64_data = split(updated_mask_data, ",")[2]
                decoded_img_bytes = base64decode(base64_data)
                
                # Load image from bytes, convert to binary, and re-save
                img_from_canvas = load(IOBuffer(decoded_img_bytes))
                binary_mask = binarize(ensure_grayscale(img_from_canvas), Otsu())

                mask_path = joinpath("public", lstrip(replace(smoothed_mask_path, r"\?.*" => ""), '/'))
                save(mask_path, binary_mask)

                smoothed_mask_path = get_timestamped_path(smoothed_mask_path)
                @info "Mask updated from canvas data and binarized."
            catch e
                @error "Failed to update mask from canvas data" exception=(e, catch_backtrace())
            end
        end
    end

    @onchange selected_folder_main begin
        if !isempty(selected_folder_main)
            is_browsing_slices = true
            is_editing_mask = false
            folder_path = joinpath("public", selected_folder_main)
            if !isdir(folder_path) 
                imgInt = ""
                msgimg = "Folder not found."
                show_verification_plot = false
                return
            end

            msi_bmp = sort(filter(f -> startswith(f, "MSI_") && endswith(f, ".bmp"), readdir(folder_path)), lt=natural)
            
            if !isempty(msi_bmp)
                current_msi = first(msi_bmp)
                imgInt = "/$(selected_folder_main)/$(current_msi)"
                text_nmass = replace(current_msi, r"MSI_|.bmp" => "")
                msgimg = "<i>m/z</i>: $(replace(text_nmass, "_" => "."))"
                
                # Use centralized plot update
                plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot(imgInt, smoothed_mask_path, imgTrans)
            else
                imgInt = ""
                msgimg = "No MSI images found in this dataset."
                show_verification_plot = false
            end
        end
    end
    
    @onbutton btn_img_plus begin
        if isempty(selected_folder_main) || is_editing_mask return end
        folder_path = joinpath("public", selected_folder_main)
        if !isdir(folder_path) return end
        msi_bmp = sort(filter(f -> startswith(f, "MSI_") && endswith(f, ".bmp"), readdir(folder_path)), lt=natural)
        new_msi = increment_image(current_msi, msi_bmp)
        if new_msi !== nothing
            current_msi = new_msi
            imgInt = "/$(selected_folder_main)/$(current_msi)"
            text_nmass = replace(current_msi, r"MSI_|.bmp" => "")
            msgimg = "<i>m/z</i>: $(replace(text_nmass, "_" => "."))"

            # Use centralized plot update
            plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot(imgInt, smoothed_mask_path, imgTrans)
        end
    end

    @onbutton btn_img_minus begin
        if isempty(selected_folder_main) || is_editing_mask return end
        folder_path = joinpath("public", selected_folder_main)
        if !isdir(folder_path) return end
        msi_bmp = sort(filter(f -> startswith(f, "MSI_") && endswith(f, ".bmp"), readdir(folder_path)), lt=natural)
        new_msi = decrement_image(current_msi, msi_bmp)
        if new_msi !== nothing
            current_msi = new_msi
            imgInt = "/$(selected_folder_main)/$(current_msi)"
            text_nmass = replace(current_msi, r"MSI_|.bmp" => "")
            msgimg = "<i>m/z</i>: $(replace(text_nmass, "_" => "."))"

            # Use centralized plot update
            plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot(imgInt, smoothed_mask_path, imgTrans)
        end
    end

    @onbutton btn_use_slice_as_mask begin
        if !isempty(imgInt)
            is_browsing_slices = false
            is_editing_mask = true
            mask_input_path = imgInt
            mask_editor_message = "Using current slice as mask input..."
            mask_editor_warning = false
            
            # Trigger processing which will update the plot via mask_input_path handler
        else
            mask_editor_message = "No slice selected. Please select a dataset and slice first."
            mask_editor_warning = true
        end
    end

    @onbutton btn_upload_mask begin
        picked_path = pick_file(; filterlist="png,bmp,jpg,jpeg")
        if !isempty(picked_path)
            target_dir = joinpath("public", "css", "masks")
            mkpath(target_dir)
            uploaded_filename = "uploaded_mask_input.png"
            destination_path = joinpath(target_dir, uploaded_filename)
            cp(picked_path, destination_path; force=true)
            
            is_browsing_slices = false
            is_editing_mask = true
            mask_input_path = get_timestamped_path("/css/masks/$uploaded_filename")
            mask_editor_message = "Uploaded mask image for processing..."
            mask_editor_warning = false
            
            # Plot will be updated via mask_input_path handler after processing
        else
            mask_editor_message = "No file selected for upload."
            mask_editor_warning = true
        end
    end

    @onbutton btn_change_slice begin
        is_browsing_slices = true
        is_editing_mask = false
        
        # Use centralized plot update
        plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot(imgInt, "", imgTrans)
        smoothed_mask_path = "" # Clear the mask path
        mask_editor_message = "Slice browsing re-enabled."
    end

    @onchange mask_input_path begin
        if isempty(mask_input_path) return end

        result = alter_image(mask_input_path, otsu_scale, noise_size_percent, hole_size_percent, smoothing_level)
        
        if result.success
            local_smoothed_path = joinpath("public", lstrip(result.path, '/'))
            slice_path_cleaned = replace(imgInt, r"\?.*" => "")
            slice_full_path = joinpath("public", lstrip(slice_path_cleaned, '/'))

            if isfile(slice_full_path) && isfile(local_smoothed_path)
                slice_img = load(slice_full_path)
                mask_img = load(local_smoothed_path)

                if size(slice_img) != size(mask_img)
                    @info "Resizing mask to match slice dimensions."
                    resized_mask = imresize(mask_img, size(slice_img))
                    save(local_smoothed_path, resized_mask)
                end
                
                smoothed_mask_path = get_timestamped_path(result.path)
                
                # Use centralized plot update
                plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot(imgInt, smoothed_mask_path, imgTrans)
                mask_editor_message = result.message
                mask_editor_warning = false
            elseif isfile(local_smoothed_path)
                smoothed_mask_path = get_timestamped_path(result.path)
                # Use centralized plot update
                plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot("", smoothed_mask_path, imgTrans)
                mask_editor_message = "Displaying generated mask. Select a dataset to overlay a slice."
                mask_editor_warning = false
            else
                mask_editor_message = "Slice or mask image not found for verification plot."
                mask_editor_warning = true
            end
        else
            mask_editor_message = result.message
            mask_editor_warning = true
        end
    end

    @onchange otsu_scale, noise_size_percent, hole_size_percent, smoothing_level begin
        if is_editing_mask
            if isempty(mask_input_path) return end

            result = alter_image(mask_input_path, otsu_scale, noise_size_percent, hole_size_percent, smoothing_level)
            
            if result.success
                local_smoothed_path = joinpath("public", lstrip(result.path, '/'))
                slice_path_cleaned = replace(imgInt, r"\?.*" => "")
                slice_full_path = joinpath("public", lstrip(slice_path_cleaned, '/'))

                if isfile(slice_full_path) && isfile(local_smoothed_path)
                    slice_img = load(slice_full_path)
                    mask_img = load(local_smoothed_path)

                    if size(slice_img) != size(mask_img)
                        @info "Resizing mask to match slice dimensions."
                        resized_mask = imresize(mask_img, size(slice_img))
                        save(local_smoothed_path, resized_mask)
                    end
                    
                    smoothed_mask_path = get_timestamped_path(result.path)
                    
                    # Use centralized plot update
                    plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot(imgInt, smoothed_mask_path, imgTrans)
                    mask_editor_message = "Mask updated."
                    mask_editor_warning = false
                elseif isfile(local_smoothed_path)
                    smoothed_mask_path = get_timestamped_path(result.path)
                    # Use centralized plot update
                    plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot("", smoothed_mask_path, imgTrans)
                    mask_editor_message = "Displaying generated mask. Select a dataset to overlay a slice."
                    mask_editor_warning = false
                else
                    mask_editor_message = "Slice or mask image not found for verification plot."
                    mask_editor_warning = true
                end
            else
                mask_editor_message = result.message
                mask_editor_warning = true
            end
        end
    end

    @onchange imgTrans begin
        if is_editing_mask
            # Use centralized plot update
            plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot(imgInt, smoothed_mask_path, imgTrans)
        end
    end

    @onbutton btn_edit_manually begin
        if !isempty(smoothed_mask_path)
            mask_path = joinpath("public", lstrip(replace(smoothed_mask_path, r"\?.*" => ""), '/'))
            if isfile(mask_path)
                mask_img = load(mask_path)
                h, w = size(mask_img)
                canvas_width = w
                canvas_height = h
                show_editor = true
            end
        else
            mask_editor_message = "No mask available for editing. Process a mask first."
            mask_editor_warning = true
        end
    end

    @onchange show_editor begin
        if !show_editor && is_editing_mask
            # When the editor dialog is closed, refresh the main plot using centralized function
            plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot(imgInt, smoothed_mask_path, imgTrans)
            is_drawing = false # Reset drawing state
        end
    end

    @onchange mouse_x, mouse_y, is_drawing begin
        if show_editor && is_drawing && mouse_x > 0 && mouse_y > 0
            try
                mask_path = joinpath("public", lstrip(replace(smoothed_mask_path, r"\?.*" => ""), '/'))
                if !isfile(mask_path) return end
                
                img = load(mask_path)
                
                actual_x = round(Int, mouse_x)
                actual_y = round(Int, mouse_y)

                if current_tool == "brush" || current_tool == "eraser"
                    if brush_size <= 0 return end
                    radius = brush_size / 2 # Treat as diameter
                    
                    for i in max(1, floor(Int, actual_y-radius)):min(size(img, 1), ceil(Int, actual_y+radius))
                        for j in max(1, floor(Int, actual_x-radius)):min(size(img, 2), ceil(Int, actual_x+radius))
                            if (i - actual_y)^2 + (j - actual_x)^2 <= radius^2
                                color = (current_tool == "brush") ? RGB(1, 1, 1) : RGB(0, 0, 0)
                                img[i, j] = color
                            end
                        end
                    end
                end
                
                save(mask_path, img)
                plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot(imgInt, smoothed_mask_path, imgTrans)
                smoothed_mask_path = get_timestamped_path(smoothed_mask_path)
                
            catch e
                @error "Editor action failed" exception=(e, catch_backtrace())
            end
        end
    end

    @onchange bucket_fill_trigger begin
        if show_editor && !isempty(bucket_fill_trigger)
            try
                mask_path = joinpath("public", lstrip(replace(smoothed_mask_path, r"\?.*" => ""), '/'))
                if !isfile(mask_path) return end
                
                img = load_and_binarize_mask(mask_path)
                
                actual_x = round(Int, bucket_fill_trigger["x"])
                actual_y = round(Int, bucket_fill_trigger["y"])

                if 1 <= actual_y <= size(img, 1) && 1 <= actual_x <= size(img, 2)
                    fill_color = Gray(1) # White
                    flood_fill!(img, actual_y, actual_x, fill_color)
                    
                    save(mask_path, img)
                    plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot(imgInt, smoothed_mask_path, imgTrans)
                    smoothed_mask_path = get_timestamped_path(smoothed_mask_path)
                    canvas_refresh_trigger = canvas_refresh_trigger[] + 1
                end
            catch e
                @error "Bucket fill failed" exception=(e, catch_backtrace())
            end
        end
    end
    
    @onchange rotate_degrees begin
        if show_editor && rotate_degrees != 0
            try
                mask_path = joinpath("public", lstrip(replace(smoothed_mask_path, r"\?.*" => ""), '/'))
                if !isfile(mask_path) return end
                
                img = load_and_binarize_mask(mask_path)
                original_size = size(img)
                
                # Normalize angle to 0, 90, 180, 270
                angle = mod(round(Int, rotate_degrees), 360)
                if angle == 0; rotate_degrees = 0; return; end

                rotated_img = if angle == 90
                    rot_r90(img)
                elseif angle == 180
                    rot180(img)
                elseif angle == 270
                    rot_l90(img)
                else
                    img # No change for other angles
                end

                if rotated_img !== img
                    # Create a new image with the original dimensions, filled with black
                    new_img = similar(img, original_size)
                    fill!(new_img, eltype(img)(0))

                    # Calculate padding/offset to center the rotated image
                    rotated_size = size(rotated_img)
                    offset_h = (original_size[1] - rotated_size[1]) ÷ 2
                    offset_w = (original_size[2] - rotated_size[2]) ÷ 2

                    # Define the region in the new image where the rotated image will be placed
                    dest_region_h = (1:rotated_size[1]) .+ offset_h
                    dest_region_w = (1:rotated_size[2]) .+ offset_w

                    # Ensure the destination region is within the bounds of the new image
                    clamped_dest_h = max(1, dest_region_h.start):min(original_size[1], dest_region_h.stop)
                    clamped_dest_w = max(1, dest_region_w.start):min(original_size[2], dest_region_w.stop)
                    
                    # Define the source region from the rotated image
                    src_h = (1:length(clamped_dest_h))
                    src_w = (1:length(clamped_dest_w))

                    if !isempty(clamped_dest_h) && !isempty(clamped_dest_w)
                        view(new_img, clamped_dest_h, clamped_dest_w) .= view(rotated_img, src_h, src_w)
                    end
                    
                    save(mask_path, new_img)
                    plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot(imgInt, smoothed_mask_path, imgTrans)
                end

                # Reset the slider and update the path
                smoothed_mask_path = get_timestamped_path(smoothed_mask_path)
                rotate_degrees = 0
                canvas_refresh_trigger = canvas_refresh_trigger[] + 1
                
            catch e
                @error "Rotation failed" exception=(e, catch_backtrace())
            end
        end
    end
    
    @onbutton btn_flip_mask begin
        if show_editor
            try
                mask_path = joinpath("public", lstrip(replace(smoothed_mask_path, r"\?.*" => ""), '/'))
                if !isfile(mask_path) return end
                
                img = load_and_binarize_mask(mask_path)
                flipped = (flip_direction == "horizontal") ? reverse(img, dims=2) : reverse(img, dims=1)
                save(mask_path, flipped)
                plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot(imgInt, smoothed_mask_path, imgTrans)
                
                smoothed_mask_path = get_timestamped_path(smoothed_mask_path)
                canvas_refresh_trigger = canvas_refresh_trigger[] + 1
                
            catch e
                @error "Flip failed" exception=(e, catch_backtrace())
            end
        end
    end

    @onchange move_mask_payload begin
        if show_editor && move_pixels > 0 && !isempty(move_mask_payload)
            try
                mask_path = joinpath("public", lstrip(replace(smoothed_mask_path, r"\?.*" => ""), '/'))
                if !isfile(mask_path) return end

                img = load_and_binarize_mask(mask_path)
                h, w = size(img)
                new_img = similar(img)
                fill!(new_img, eltype(img)(0))

                dx, dy = 0, 0
                direction = move_mask_payload["direction"]
                if direction == "right"; dx = move_pixels; end
                if direction == "left"; dx = -move_pixels; end
                if direction == "down"; dy = move_pixels; end
                if direction == "up"; dy = -move_pixels; end

                src_y_range = max(1, 1-dy):min(h, h-dy)
                src_x_range = max(1, 1-dx):min(w, w-dx)
                dest_y_range = max(1, 1+dy):min(h, h+dy)
                dest_x_range = max(1, 1+dx):min(w, w+dx)
                
                if !isempty(src_y_range) && !isempty(src_x_range) && !isempty(dest_y_range) && !isempty(dest_x_range)
                    view(new_img, dest_y_range, dest_x_range) .= view(img, src_y_range, src_x_range)
                end

                save(mask_path, new_img)
                plotdata_verify, plotlayout_verify, show_verification_plot = update_main_plot(imgInt, smoothed_mask_path, imgTrans)
                smoothed_mask_path = get_timestamped_path(smoothed_mask_path)
                canvas_refresh_trigger = canvas_refresh_trigger[] + 1

            catch e
                @error "Move failed" exception=(e, catch_backtrace())
            end
        end
    end
    
    @onbutton btn_save_final_mask begin
        try
            parent_folder = split(selected_folder_main, '/')[end]
            final_mask_name = "$(parent_folder).png"
            
            source_path = joinpath("public", lstrip(replace(smoothed_mask_path, r"\?.*" => ""), '/'))
            target_dir = joinpath("public", "css", "masks")
            mkpath(target_dir)
            final_path = joinpath(target_dir, final_mask_name)
            
            if isfile(source_path)
                # Copy the mask to final location
                cp(source_path, final_path; force=true)
                
                # Update registry directly in the handler
                reg_path = abspath(joinpath(@__DIR__, "public", "registry.json"))
                registry = isfile(reg_path) ? JSON.parsefile(reg_path) : Dict{String, Any}()

                full_mask_path = abspath(final_path)
                
                # Update the registry entry
                if haskey(registry, selected_folder_main)
                    registry[selected_folder_main]["mask_path"] = full_mask_path
                    registry[selected_folder_main]["has_mask"] = true
                else
                    # Create a new entry if folder doesn't exist in registry
                    registry[selected_folder_main] = Dict(
                        "mask_path" => full_mask_path,
                        "has_mask" => true,
                        "is_imzML" => true,
                        "processed_date" => string(Dates.now())
                    )
                end
                
                # Save the updated registry
                open(reg_path, "w") do f
                    JSON.print(f, registry, 4)
                end
                
                @info "Registry updated with mask: $(final_mask_name)"
                
                # Update UI state
                mask_editor_message = "Final mask saved as $(final_mask_name)"
                mask_editor_warning = false
                show_editor = false
            else
                mask_editor_message = "Source mask file not found: $(source_path)"
                mask_editor_warning = true
            end
            
        catch e
            @error "Save final mask failed" exception=(e, catch_backtrace())
            mask_editor_message = "Error saving final mask: $(sprint(showerror, e))"
            mask_editor_warning = true
        end
    end

    @onchange isready begin
        if isready && !registry_init_done
            sleep(1.0) # Give frontend time to initialize
            try
                println("Synchronizing registry for mask editor...")
                reg_path = abspath(joinpath(@__DIR__, "public", "registry.json"))
                # Assuming load_registry is available from MSI_src or a similar utility file
                # For now, handle its absence gracefully if it's not explicitly defined here.
                registry = isfile(reg_path) ? load_registry(reg_path) : Dict{String, Any}()

                public_dirs = isdir("public") ? readdir("public") : []
                ignored_dirs = ["css", "masks"]
                
                dataset_dirs = filter(d -> isdir(joinpath("public", d)) && !(d in ignored_dirs), public_dirs)
                
                registry_keys = Set(keys(registry))
                folder_set = Set(dataset_dirs)

                new_folders = setdiff(folder_set, registry_keys)
                for folder in new_folders
                    println("Found new folder: $folder")
                    registry[folder] = Dict(
                        "source_path" => "unknown (manually added)",
                        "processed_date" => "unknown",
                        "metadata" => Dict(),
                        "is_imzML" => true # Assume folder contains images if found this way
                    )
                end

                removed_folders = setdiff(registry_keys, folder_set)
                for folder in removed_folders
                    delete!(registry, folder)
                end

                if !isempty(new_folders) || !isempty(removed_folders)
                    println("Registry changed, saving...")
                    open(reg_path, "w") do f
                        JSON.print(f, registry, 4)
                    end
                end
                
                all_folders = sort(collect(keys(registry)), lt=natural)
                img_folders = filter(folder -> get(get(registry, folder, Dict()), "is_imzML", false), all_folders)

                available_folders = deepcopy(all_folders)
                image_available_folders = deepcopy(img_folders)
                println("Mask editor UI lists updated. All: $(length(available_folders)), Images: $(length(image_available_folders))")
            catch e
                @warn "Mask editor registry synchronization failed: $e"
                available_folders = []
                image_available_folders = []
            finally
                registry_init_done = true
            end
        end
    end
    
    @methods """
    // Close the default methods object and define our own component structure
    },

    created() {
        // Initialize non-reactive data properties here
        this.canvas = null;
        this.ctx = null;
        this.isDrawingOnCanvas = false;
        this.lastX = 0;
        this.lastY = 0;
        this.imgObj = null;
    },

    watch: {
        show_editor(newValue) {
            if (newValue) {
                // Wait for the dialog to render before initializing canvas
                setTimeout(() => { this.initCanvas() }, 100);
            }
        },
        canvas_refresh_trigger(newValue) {
            if (newValue > 0) {
                this.initCanvas();
            }
        }
    },

    methods: {
        // Re-opened methods object for all our functions
        initCanvas() {
            this.canvas = document.getElementById('maskCanvas');
            if (!this.canvas) {
                console.error("Canvas element not found!");
                return;
            }
            this.ctx = this.canvas.getContext('2d');

            this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

            this.imgObj = new Image();
            this.imgObj.onload = () => {
                this.ctx.drawImage(this.imgObj, 0, 0, this.canvas.width, this.canvas.height);
            };
            this.imgObj.src = this.smoothed_mask_path.split('?')[0] + '?t=' + Date.now();
        },

        startDrawing(event) {
            if (this.current_tool === 'drag' || !this.ctx) return;

            this.isDrawingOnCanvas = true;
            const { x, y } = this.getCanvasMousePosition(event);
            this.lastX = x;
            this.lastY = y;

            if (this.current_tool === 'bucket') {
                this.bucket_fill_trigger = { x: x, y: y, t: Date.now() };
                this.isDrawingOnCanvas = false;
                return;
            }
            
            this.ctx.beginPath();
            this.ctx.moveTo(this.lastX, this.lastY);
        },

        draw(event) {
            if (!this.isDrawingOnCanvas || this.current_tool === 'drag' || this.current_tool === 'bucket' || !this.ctx) return;

            const { x, y } = this.getCanvasMousePosition(event);

            this.ctx.lineWidth = this.brush_size;
            this.ctx.lineCap = 'round';
            this.ctx.lineJoin = 'round';

            if (this.current_tool === 'brush') {
                this.ctx.strokeStyle = 'white';
                this.ctx.globalCompositeOperation = 'source-over';
            } else if (this.current_tool === 'eraser') {
                this.ctx.strokeStyle = 'black';
                this.ctx.globalCompositeOperation = 'source-over';
            }

            this.ctx.lineTo(x, y);
            this.ctx.stroke();
        },

        stopDrawing() {
            if (this.isDrawingOnCanvas) {
                this.isDrawingOnCanvas = false;
                this.ctx.closePath();
                this.syncMaskToServer();
            }
        },

        getCanvasMousePosition(event) {
            const rect = this.canvas.getBoundingClientRect();
            const scaleX = this.canvas.width / rect.width;
            const scaleY = this.canvas.height / rect.height;
            const x = (event.clientX - rect.left) * scaleX;
            const y = (event.clientY - rect.top) * scaleY;
            return { x, y };
        },

        syncMaskToServer() {
            if (!this.canvas) return;
            const dataURL = this.canvas.toDataURL('image/png');
            this.updated_mask_data = dataURL;
        },

        // Kept for compatibility with bucket tool logic which uses original image coordinates
        updateMousePosition(event) {
            const rect = event.target.getBoundingClientRect();
            const scale = this.editor_scale || 1;
            this.mouse_x = (event.clientX - rect.left) / scale;
            this.mouse_y = (event.clientY - rect.top) / scale;
        }
    // The closing brace for methods is intentionally omitted, as Genie/Stipple will add it.
    """
    
end

@page("/mask", "mask.jl.html")

end # module
