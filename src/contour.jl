# Creates contour lines overlay from elevation data.
# Output is image files for manual touch-up.

"""
    contour_lines_overlay(sb::SheetBuilder)
    contour_lines_overlay(fofo)

"""
contour_lines_overlay(sb::SheetBuilder) = contour_lines_overlay(full_folder_path(sb))
function contour_lines_overlay(fofo)
    if isfile(joinpath(fofo, CONTOUR_FNAM))
        @debug "$CONTOUR_FNAM in $fofo already exists. Exiting `contour_lines_overlay`."
        return true
    end
    _contour_lines(fofo)
end
function _contour_lines(fofo)
    true
end