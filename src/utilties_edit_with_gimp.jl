# Utilty functions
#
# Functions for
# - aiding manual touch-up of elevation .tif files
# - opening other image files or image matrices for editing or inspection
#
# Gimp is called as an external process or task.


"""
    edit_in_gimp(sb::SheetBuilder, mini, maxi)
    ---> Matrix

Use case:

If you notice, typically, errors in the data along the edges of a lake that the pipeline couldn't automatically fix, manually
editing the elevations close to that lake can be useful.

Manually modify cells in the elevation file corresponding to `sb`, and with elevations in the interval [`mini`, `maxi`>.

WARNING: This could easily result in unwanted degradation, so take note that:

If the elevation matrix contains Float32 values, it means that it distinguishes between elevation
12.345678 and 12.345679, although the measurements are not that precise. Since the `BitmapMaps` pipeline uses
curvature and other tricks to detect terrain features or glitches in the raw data, changing the numeric
resolution may have consequences and should be avoided.

Although `gimp` has some ability to deal with grayscale Float32, it's practically limited
to easily working with grayscale numbers 0 (black) to 255 (white).

DO:
1) In the gimp editor, delete any pixels in areas that you're not interested in changing. Deleted pixels will not
   affect changes in the elevation file.
2) Don't use unnecessarily large intervals `mini` to `maxi`. Five metres is usually enough to easily fix glitches.


# Example

This opens the consolidated elevation file for the SheetBuilder. Every pixel not in [1, 2> will be transparent and zero-valued.

When user saves, then closes the editor, every pixel not zero-valued will overwrite the corresponding elevation in the consolidated elevation file.
```
julia> edit_in_gimp(sb, 1f0, 2f0)
Reading values from Consolidated.tif
 Warning, this is based on tifdic.jls only
Folder for consolidated file: C:\\Users\\f\\BitmapMaps/1 1 18145 6894028 71625 6971028\\7-7__63985-6960028__71625-6971028
Consolidated file internal bounding box, utm coordinates: (min_x = 63985, min_y = 6960030, max_x = 71623, max_y = 6971028)
Probable source files for this bounding box: (min_x = 63985, min_y = 6960030, max_x = 71623, max_y = 6971028)
C:\\Users\\f\\BitmapMaps\\1 1 18145 6894028 71625 6971028\\33-111-135.tif    Bounding box (min_x = 65425, min_y = 6950995, max_x = 80435, max_y = 6966005))
C:\\Users\\f\\BitmapMaps\\1 1 18145 6894028 71625 6971028\\33-111-136.tif    Bounding box (min_x = 65425, min_y = 6965995, max_x = 80435, max_y = 6981005))
C:\\Users\\f\\BitmapMaps\\1 1 18145 6894028 71625 6971028\\33-110-135.tif    Bounding box (min_x = 50425, min_y = 6950995, max_x = 65435, max_y = 6966005))
C:\\Users\\f\\BitmapMaps\\1 1 18145 6894028 71625 6971028\\33-110-136.tif    Bounding box (min_x = 50425, min_y = 6965995, max_x = 65435, max_y = 6981005))
Encoding to GrayA.
1.0 m =   0.0 %        0 / 255
1.2 m =  25.0 %       64 / 255
1.5 m =  50.0 %      128 / 255
1.8 m =  75.0 %      191 / 255
2.0 m =   0.0 %        0 / 255
Saving scratch file C:\\Users\\F\\AppData\\Local\\Temp\\jl_M3wRgjsf2N.png
Opening scratch file in gimp. This process is blocking until Gimp is closed.
 If other files are already open, returns immediately
 ERROR: "No edit change registered - in Gimp, press Alt+F, W, Ctrl+Q, Ctrl+D to re-export, then exit. C:\\\\Users\\\\F\\\\AppData\\\\Local\\\\Temp\\\\jl_M3wRgjsf2N.png"
```
"""
function edit_in_gimp(sb::SheetBuilder, mini, maxi)
    printstyled("Reading values from $CONSOLIDATED_FNAM \n ", color = :green)
    z = elevation_full_res(sb; display_sources = true)
    edit_in_gimp!(z, mini, maxi)
    modify_consolidated_file(sb, permutedims(z))
    z
end

"""
    edit_in_gimp!(M::Matrix{S}, mini::S, maxi::S) where S
    ---> typeof(M)

Modifes M in-place based on user edits made in the editor `gimp`.
Only non-zero-valued pixel will affect changes in M, so deleting pixels in the external
editor is harmless.
"""
function edit_in_gimp!(M::Matrix{S}, mini::S, maxi::S) where S
    printstyled("Encoding to GrayA. \n", color = :green)
    img = encode_GrayA(M, mini, maxi)
    ffnam = tempname() * ".png"
    printstyled("Saving scratch file $ffnam \n", color = :green)
    save_png_with_phys(ffnam, img)
    # This waits for the external process to end. Use ctrl-c to interrupt.
    # Note that if gimp is already open, the process will return at once.
    # If not, it waits for gimp.exe to close.
    printstyled("Opening scratch file in gimp. This process is blocking until Gimp is closed. \n If other files are already open, returns immediately\n ", color = :green)
    open_in_gimp(ffnam)
    # Let's extract from the green channel (red and blue can be changed or not, no matter)
    mod = load(ffnam)
    if mod == img
        throw("No edit change registered - in Gimp, press Alt+F, W, Ctrl+Q, Ctrl+D to re-export, then exit. $ffnam")
    else
        printstyled("Decoding GrayA to $S. The α channel is ignored.\n ", color = :green)
        decode_GrayA!(M, mod, mini, maxi)
    end
    M
end


"""
    modify_consolidated_file(sb, M)
    --> String


Provided that the SheetBuilder sb has a ready made consolidated geoarray file,
and that it has the same size and shape as M, replace the cell values in that
file to match M. The old consolidated file is written over.

Returns the file name of the consolidated file.
"""
function modify_consolidated_file(sb, M)
    ffnam = joinpath(full_folder_path(sb), CONSOLIDATED_FNAM)
    @assert isfile(ffnam) ffnam
    modify_geoarray_file(ffnam, M)
end



"""
    open_as_temp_in_gimp(img)
    ---> Task

Save the image `img` to a temporary file, and opens that file in `gimp` editor.

Non-blocking. If `gimp` was already running, returns a done task. If not, returns
a started task, which will be finished once `gimp` is closed.
"""
function open_as_temp_in_gimp(img)
    @async let
        fnam = tempname()
        save_png_with_phys(fnam, img)
        open_in_gimp(fnam)
    end
end


"""
    open_in_gimp(fnam)
    ---> Process

If `gimp` is not already open: Opens file name 'fnam' in the image editor `gimp` in a blocking manner.

If already open, possibly with another file, the process will exit immediately, and `gimp` will continue to run.
"""
function open_in_gimp(fnam)
    # TODO Consider: Make the path to the application configureable.
    #     If so, also change the function name.
    # Input check
    isfile(fnam) || throw(ErrorException("File does not exists: $fnam"))
    # Determine the operating system
    if Sys.iswindows()
        # Common GIMP installation paths on Windows
        gimp_paths = [
            joinpath("C:\\", "Program Files", "GIMP 2", "bin", "gimp-2.10.exe"),
            joinpath("C:\\", "Program Files", "GIMP 2", "bin", "gimp-2.8.exe"),
            joinpath("C:\\", "Program Files", "GIMP 2", "bin", "gimp.exe"),
        ]
    elseif Sys.islinux()
        # Common GIMP installation paths on Linux
        gimp_paths = [
            "/usr/bin/gimp",
            "/usr/local/bin/gimp",
        ]
    elseif Sys.isapple()
        # Common GIMP installation paths on macOS
        gimp_paths = [
            "/Applications/GIMP.app/Contents/MacOS/gimp",
            "/Applications/GIMP-2.10.app/Contents/MacOS/gimp",
        ]
    else
        error("Unsupported operating system")
    end
    # Try to find the GIMP executable
    gimp_executable = nothing
    for path in gimp_paths
        if isfile(path)
            gimp_executable = path
            break
        end
    end
    if gimp_executable === nothing
        error("GIMP executable not found. Please ensure GIMP is installed.")
    end
    # Run GIMP with the specified file
    run(`$gimp_executable $fnam`)
end






"""
    modify_geoarray_file(ffnam, M::Matrix{T}) where T
    --> String

Provided that `ffnam` is the file name of a geoarray of the same size and shape as M,
replace the cell values in that file to match M. The old file is written over.

Returns ffnam.
"""
function modify_geoarray_file(ffnam, M::Matrix{T}) where T
    g = readclose(ffnam)
    modify_geoarray!(g, M)
    # Write modified file
    printstyled("Writing modified geoarrays file $ffnam \n", color = :green)
    GeoArrays.write(ffnam, g)
end


"""
    modify_geoarray!(g::GeoArray{T, 2, Matrix{T}}, M::Matrix{T}) where T <: Float32
    ---> Vector{Float64}

"""
function modify_geoarray!(g::GeoArray{T, 2, Matrix{T}}, M::Matrix{T}) where T <: Float32
    @assert size(g) == size(M) "Forgot permutedims? $(size(g)[1:2]) !==  $(size(M))"
    g.A[:, :] .= M
end



"""
    scaleminmax00(min::T, max::T) where T

Like `ImageCore.scaleminmax`, return a function f which maps values less than or equal to min to 0,
linearly increasing towards 1 at max. At or above max, also returns 0.
"""
function scaleminmax00(min::T, max::T) where T
    @inline function(x)
        xp, minp, maxp = promote(x, min, max)  # improves performance to promote explicitly
        y = clamp(xp, minp, maxp)
        (y == maxp || y == minp) ? zero(y) : (y-minp)/(maxp-minp)
    end
end


"""
    encode_GrayA(M::Matrix{S}, mini, maxi) where S <: Union{Float64, Float32, Bool, Int64}

Takes potentially high-resolution data in matrix form. If values fall in the interval [mini, maxi>,
will scale those linearly to 255 individual gray tones. Otherwise, returns zero.

Also adds an alpha (transparency channel), which takes values 1 where there is data in the interval, and 0
elsewhere.

If max and min are far apart, will reduce the value resolution accordingly.
"""
function encode_GrayA(M::Matrix{S}, mini, maxi) where S <: Union{Float64, Float32, Bool, Int64}
    if maxi - mini > 60
        @warn "Large difference between max and min may reduce the vertical resolution too much, the resolution is (maxi-mini) /256 !"
    end
    f = scaleminmax00(S(mini), S(maxi))
    rng = range(mini, maxi, length = 5)
    for x in rng
        printstyled(round(x, digits = 1), " m = ", lpad(round(f(x) * 100), 5), " %    ", lpad(Int(round(f(x) * 255)), 5), " / 255\n", color = :green)
    end
    map(z -> GrayA{N0f8}(f(z), f(z) > 0), M)
end

"""
    decode_GrayA!(M::Matrix{S}, img::Matrix{GrayA{T}}, mini::S, maxi::S) where {S, T}

Modify in-place values from M, provided the equally sized matrix img has a correspondingly placed non-zero value.

Maps 0 to `mini`, linearly up to 1 to `maxi`.

The α channel is ignored.
"""
function decode_GrayA!(M::Matrix{S}, img::Matrix{GrayA{T}}, mini::S, maxi::S) where {S, T}
    @assert size(M) == size(img)
    for I in CartesianIndices(axes(M))
        c = gray(img[I]) * (maxi - mini) + mini
        if c !== mini
            M[I] = c
        end
    end
    M
end


