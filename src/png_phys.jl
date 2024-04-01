
"""
    save_png_with_phys(fpath::String,
    image::S,
    pt_per_m::Int;
    compression_level::Integer = Z_BEST_SPEED,
    compression_strategy::Integer = Z_RLE,
    filters::Integer = Int(PNG_FILTER_PAETH),
    file_gamma::Union{Nothing,Float64} = nothing,
    background::Union{Nothing,UInt8,AbstractGray,AbstractRGB} = nothing,
) where {
    T,
    S<:Union{AbstractMatrix{T},AbstractArray{T,3}}
}
    ---> nothing

Set the pHYs chunk values to pt_per_m for both vertical and horizontal, unlike 'PNGFiles.save'.
This is used for printing with e.g. IrfanView, but is ignored by the default win64 printing
via default apps.
"""
function save_png_with_phys(fpath::String,
        image::S,
        density_pt_per_m::Int;
        compression_level::Integer = Z_BEST_SPEED,
        compression_strategy::Integer = Z_RLE,
        filters::Integer = Int(PNG_FILTER_PAETH),
        file_gamma::Union{Nothing,Float64} = nothing,
        background::Union{Nothing,UInt8,AbstractGray,AbstractRGB} = nothing,
        ) where {
            T,
            S<:Union{AbstractMatrix{T},AbstractArray{T,3}}
        }
    @assert Z_DEFAULT_STRATEGY <= compression_strategy <= Z_FIXED
    @assert Z_NO_COMPRESSION <= compression_level <= Z_BEST_COMPRESSION
    @assert 2 <= ndims(image) <= 3
    @assert size(image, 3) <= 4

    fp = ccall(:fopen, Ptr{Cvoid}, (Cstring, Cstring), fpath, "wb")
    fp == C_NULL && error("Could not open $(fpath) for writing")

    png_ptr = create_write_struct()
    @debug "save_png_with_phys:" fpath png_ptr pt_per_m
    info_ptr = create_info_struct(png_ptr)
    png_init_io(png_ptr, fp)

    # Set the pHYs chunk here, unlike 'PNGFiles.save'
    res_x = png_uint_32(density_pt_per_m)
    res_y = png_uint_32(density_pt_per_m)
    unit_type = Cint(1)
    
    png_set_pHYs(png_ptr, info_ptr, res_x, res_y, unit_type)
    
    # Continue as with normal save
    _save(png_ptr, info_ptr, image,
        compression_level=compression_level,
        compression_strategy=compression_strategy,
        filters=filters,
        file_gamma=file_gamma,
        background=background
    )
    close_png(fp)
    return
end

"""
    get_pHYs_chunk_res_x_y_unit(fpath; silent = false)
    ---> res_x::Int, res_y, unit_type

unit_type == 1 signifies res_x and res_y are given in units
pixels per printed meter.

For inspection / checking
"""
function get_pHYs_chunk_res_x_y_unit(fpath; silent = false)
    fp = open_png(fpath) # pointer
    png_ptr = create_read_struct()
    info_ptr = create_info_struct(png_ptr)
    png_init_io(png_ptr, fp)
    # Since we don't really know much of what we're doing, we're keeping code from _inspect_png_read
    # which just might make the inspection more generic.
    png_set_sig_bytes(png_ptr, PNG_BYTES_TO_CHECK)
    png_read_info(png_ptr, info_ptr)

    width = png_get_image_width(png_ptr, info_ptr)
    height = png_get_image_height(png_ptr, info_ptr)

    color_type_orig = png_get_color_type(png_ptr, info_ptr)
    color_type = color_type_orig
    bit_depth_orig = png_get_bit_depth(png_ptr, info_ptr)
    bit_depth = bit_depth_orig
    backgroundp = png_color_16p()
    if png_get_bKGD(png_ptr, info_ptr, Ref(backgroundp)[]) != 0
        png_set_background(png_ptr, backgroundp, PNG_BACKGROUND_GAMMA_FILE, 1, 1.0)
    end

    if color_type == PNG_COLOR_TYPE_PALETTE
        png_set_palette_to_rgb(png_ptr)
        color_type = PNG_COLOR_TYPE_RGB
    end

    if color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8
        png_set_expand_gray_1_2_4_to_8(png_ptr)
        png_set_packing(png_ptr)
        bit_depth = UInt8(8)
    end

    # Now do the necessary part
    bit_depth == 16 && png_set_swap(png_ptr)
    res_x = Ref{png_uint_32}(0)
    res_y = Ref{png_uint_32}(0)
    unit_type = Ref{Cint}(0)
    png_get_pHYs(png_ptr, info_ptr, res_x, res_y, unit_type)
    # Human-readable feedback
    if ! silent
        if unit_type[] == 1
            w_mm = 1000 * width / res_x[]
            h_mm = 1000 * height / res_y[]
            dens_x_dpi = res_x[] * 25.4 / 1000
            dens_y_dpi = res_x[] * 25.4 / 1000
            println("Print width = $w_mm mm from png:pHYs")
            println("Print height = $h_mm mm from png:pHYs")
            println("Print resolution x = $dens_x_dpi dpi from png:pHYs")
            println("Print resolution y = $dens_y_dpi dpi from png:pHYs")
        else
            @show res_x[] res_y[] unit_type[]
        end
        println("Display width = $(Int(width))")
        println("Display height = $(Int(height))")
    end
    png_destroy_read_struct(Ref{Ptr{Cvoid}}(png_ptr), Ref{Ptr{Cvoid}}(info_ptr), C_NULL)
    close_png(fp)
    res_x[], res_y[], unit_type[] 
end