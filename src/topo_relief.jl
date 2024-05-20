# This renders a topographic relief map.
"""
    topo_relief(sb::SheetBuilder)
    topo_relief(fofo)

"""
topo_relief(sb::SheetBuilder) = topo_relief(full_folder_path(sb))
function topo_relief(fofo) # No, we'll use the iterators and so on, not just the folder.
    if isfile(joinpath(fofo, TOPORELIEF_FNAM))
        @debug "$TOPORELIEF_FNAM in $fofo already exists. Exiting `topo_relief`."
        return true
    end
    if ! isfile(joinpath(fofo, HYPSOMETRIC_FNAM))
        @debug "$HYPSOMETRIC_FNAM in $fofo does not exist. Exiting `topo_relief`."
        return false
    end
    _topo_relief(fofo)
end
function _topo_relief(fofo)
    true
end