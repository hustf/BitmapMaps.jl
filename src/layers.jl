# Merges foreground layers onto background using the α / opacity channel.
# Output is an image file per sheet.

"""
    join_layers(sb::SheetBuilder)
    join_layers(fofo)

"""
join_layers(sb::SheetBuilder) = join_layers(full_folder_path(sb))
function join_layers(fofo)
    if isfile(joinpath(fofo, CONTOUR_FNAM))
        @debug "$CONTOUR_FNAM in $fofo already exists. Exiting `join_layers`."
        return true
    end
    _join_layers(fofo)
end
function _join_layers(fofo)
    true
end





# TODO change things
function add_lakes_and_save_to_routemap(i, j)
    CONSOLIDATED_FNAM = "Consolidated.tif"
    TOPORELIEF_FNAM = "Toporelief.png"
    WATER_FNAM = "Water.png"
    CONTOUR_FNAM = "Contour.png"
    relief_fnam =  joinpath(homedir(), "Høydedata\\$i $j elev lines downsampled by factor 3.png")
    lake_fnam = joinpath(homedir(), "Høydedata\\Lakes\\$i $j lakes.png")
    dest_fnam = joinpath(homedir(), ".julia\\dev\\RouteMap\\example\\Split\\$i $j bg.png")
    relief_img = Images.load(relief_fnam)
    lake_img = Images.load(lake_fnam)
    mask = [Bool(round(alpha(p), digits = 0)) for p in lake_img]
    comp_img = pick_by_mask(lake_img, relief_img, mask)
    saveshow(dest_fnam, comp_img)
end
