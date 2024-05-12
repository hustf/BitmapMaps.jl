# Creates water overlay from elevation data.
# Output is image files for manual touch-up.

"""
    water_overlay(sb::SheetBuilder)
    water_overlay(fofo)

"""
water_overlay(sb::SheetBuilder) = water_overlay(full_folder_path(sb))
function water_overlay(fofo)
    if isfile(joinpath(fofo, WATER_FNAM))
        @debug "$WATER_FNAM in $fofo already exists. Exiting `water_overlay`."
        return true
    end
    _water_overlay(fofo)
end
function _water_overlay(fofo)
    true
end