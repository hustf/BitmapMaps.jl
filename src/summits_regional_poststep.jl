# A 'regional' or 'postprocessing' step in the pipeline,
# which requires that all `summits_on_sheet` have been run for all sheets.
# Reduces the summits list based on prominence, and harvests geographical names.
#
# The regionally calculated summit prominence and position is used to:
# 1) Update Summits.csv and Markers.png (in this step)
# 2) From the calling context(the pipeline), re-run `make_vector_graphics` and `join_layers`.
#

"""
    summits_regional_update(smb::SheetMatrixBuilder, ffna_graph)
    ---> Bool

1) Harvest summit prominence & etc. from the regional elevation graph.
2) Reduces the list of summits based on prominence.
3) Retrieves names for the shorter list of summits
4) Distributes info to Summits.csv
5) Distributes info to Markers.png
"""
function summits_regional_update(smb::SheetMatrixBuilder, ffna_graph)
    # Early exit
    isfile(ffna_graph) || throw(ErrorException("$ffna_graph is missing"))
    graphtime = mtime(ffna_graph)
    summitsfiles = [joinpath(full_folder_path(sb), SUMMITS_FNAM) for sb in smb]
    summits_newer_than_regional_graph  = mtime.(summitsfiles) .> graphtime
    if all(summits_newer_than_regional_graph  )
        @debug "    The regional elevation graph $(splitpath(ffna_graph)[end]) exists and is older than all $(SUMMITS_FNAM). Exiting `summits_regional_update`"
        return true
    else
        # Some SUMMITS_FNAM are older than the regional graph.
        # Still, if those are empty (no summits at all in those sheets), we can early exit.
        sheet_nos_of_older = findall(iszero, summits_newer_than_regional_graph )
        sheet_nos_of_older_with_content = filter(sheet_nos_of_older) do shno
            ffna_sum = summitsfiles[shno]
            summits_data = readdlm(ffna_sum, '\t')[2:end,:]
            size(summits_data)[1] !== 0
        end
        if isempty(sheet_nos_of_older_with_content)
            @debug "    The regional elevation graph $(splitpath(ffna_graph)[end]) exists and is older than all $(SUMMITS_FNAM) with content. Exiting `summits_regional_update`"
            return true
        end
    end
    # Configuration pararameters
    promlev_prom = get_config_value("Markers", "Prominence level [m], prominent summit", Int)
    promlev_obsc = get_config_value("Markers", "Prominence level [m], obscure summit", Int)
    prom_levels = [promlev_obsc, promlev_prom]
    symbol_prom = get_config_value("Markers", "Symbol prominent summit", String)
    symbol_obsc = get_config_value("Markers", "Symbol obscure summit", String)
    symbol_size_prom = get_config_value("Markers", "Size prominent summit symbol", Int)
    symbol_size_obsc = get_config_value("Markers", "Size obscure summit symbol", Int)
    summit_symbols = [symbol_obsc, symbol_prom]
    symbol_sizes = [symbol_size_obsc, symbol_size_prom]
    # Get the regional graph.
    gr = get_graph(ffna_graph)
    @assert bbox_internal(gr) == bbox_internal(smb)
    # Get elevation, prominence, position in utm from the graph.
    # A few summits will be 'doubled up' due to rounding & etc.
    # For such 'fake' summits, we will not find a corresponding 'σ' stress
    # value in the SUMMITS_FNAM file.
    @debug "    Harvest summits data from regional graph, at $(nowstring())"
    vz, vutm, vprom, vsaddle, vtaller = harvest_summits_data_from_graph(gr)
    # Get sheet indices, cell indices, names
    @debug "    Harvest other summits data, including online if config allows, at $(nowstring())"
    vsheet_ij, vcell_ij, vσ, vname  = harvest_other_summits_data(smb, vutm)
    # Place the vectors in a dictionary indexed by vsheet_ij.
    @debug "    Distribute summits data to summits text files and markers layer, at $(nowstring())"
    # We drop the vectors which were not present in the preliminary SUMMITS_FNAM
    # file. Those now have σ = NaN32.
    # We also drop the few summits which were included in the preliminary SUMMITS_FNAM
    # because their dominating summit was on a neighbouring sheet, but have too low
    # prominence when calculated from the regional elevation graph.
    dic = Dict{Tuple{Int64, Int64}, Vector{Any}}()
    for ij in unique(vsheet_ij)
        for (sheet_ij, z, prom, utm, cell_ij, σ, name) in zip(vsheet_ij, vz, vprom, vutm, vcell_ij, vσ, vname)
            if sheet_ij == ij
                if ! isnan(σ)
                    if prom >= promlev_obsc
                        if haskey(dic, ij)
                            push!(dic[ij], [z, prom, utm, cell_ij, σ, name])
                        else
                            push!(dic, ij => [[z, prom, utm, cell_ij, σ, name]])
                        end
                    end
                end
            end
        end
    end
    #
    # Distribute the collected data
    #
    distribute_summits_data(smb, dic, prom_levels, summit_symbols, symbol_sizes)
end

"""
    distribute_summits_data(smb, dic, prom_levels, summit_symbols, symbol_sizes)
"""
function distribute_summits_data(smb, dic, prom_levels, summit_symbols, symbol_sizes)
    #
    # Write top-level csv. This is for easy inspection.
    #
    # Note the name is not SUMMITS_FNAM, since we want to be able to extend the number of sheets
    # in a folder. The name instead reflects the size of the regional map.
    ffnam = joinpath(full_folder_path(smb), bbox_external_string(smb) * ".csv")
    # Unpack all sheets for top-level vectors
    v = vcat(values(dic)...) # Nested vector{Any}
    vz, vprom, vutm, vcell_ij, vσ, vname = [getindex.(v, i) for i in 1:6]
    # Find an order: First prominent summits, then other summits.
    # Otherwise, ordered by height z.
    sortval = [z + (p > prom_levels[2] ? 10000 : 0) for (z, p) in zip(vz, vprom)]
    order = sortperm(sortval; rev = true)
    # Names will be used as headers in the .csv file
    # nt :: @NamedTuple{vz::Vector{Float32}, vprom::Vector{Int64}, vutm::Vector{Tuple{Int64, Int64}}, vcell_ij::Vector{Tuple{Int64, Int64}}, vσ::Vector{AbstractFloat}, vname::Vector{String}}
    nt = (; Elevation_m = vz[order], Prominence_m = vprom[order], Utm = vutm[order], Cell_index = vcell_ij[order], Stress = vσ[order], Name = vname[order] )
    write_named_tuple_to_csv(ffnam, nt)
    #
    # Write SUMMITS_FNAM in every folder. These are the ones used for labelling the .svg.
    #
    for (sheet_ij, v) in dic
        sb = smb[sheet_ij...]
        ffnam = joinpath(full_folder_path(sb), SUMMITS_FNAM)
        # The rest of the loop is a copy of writing the top-level .csv above.
        vz, vprom, vutm, vcell_ij, vσ, vname = [getindex.(v, i) for i in 1:6]
        sortval = [z + (p > prom_levels[2] ? 10000 : 0) for (z, p) in zip(vz, vprom)]
        order = sortperm(sortval; rev = true)
        nt = (; Elevation_m = vz[order], Prominence_m = vprom[order], Utm = vutm[order], Cell_index = vcell_ij[order], Stress = vσ[order], Name = vname[order] )
        write_named_tuple_to_csv(ffnam, nt)
    end
    #
    # Update the preliminary file MARKERS_FNAM which was generated by `summits_on_sheet`.
    #
    for (sheet_ij, v) in dic
        sb = smb[sheet_ij...]
        mffna = joinpath(full_folder_path(sb), MARKERS_FNAM)
        vz, vprom, vutm, vcell_ij, vσ, vname = [getindex.(v, i) for i in 1:6]
        # Type change for convenience and speed
        vI = CartesianIndex.(vcell_ij)
        # Make an empty matrix, output image size
        prom_mat = fill(0f0, axes(sb.cell_iter));
        # Set the values at vcell_ij to the corresponding prominence values
        prom_mat[vI] .= Float32.(vprom)
        # Make an image with prominent and obscure summit symbols as specified in the configuration
        bwres = draw_summit_marks(prom_mat, vI, prom_levels, summit_symbols, symbol_sizes)
        # Feedback
        display_if_vscode(bwres)
        @debug "    Saving updated $mffna"
        # Make a transparent color image, save it.
        save_png_with_phys(mffna, map(bwres) do pix
            pix == true && return RGBA{N0f8}(0., 0, 0, 1)
            RGBA{N0f8}(0., 0, 0, 0)
        end)
    end
    true
end


function harvest_other_summits_data(smb, vutm)
    # Get indices of which sheet number each summit belongs to
    f_utm_to_sheet_index = func_utm_to_sheet_index(smb)
    vsheet_ij =  map(f_utm_to_sheet_index, vutm) # Tuples
    vsb = [smb[sij...] for (sij, utm) in zip(vsheet_ij, vutm)] # Sheet builders
    # Get indices of which cell index each summit belongs to
    vcell_ij = [func_utm_to_cell_index(sb)(utm) for (sb, utm) in zip(vsb, vutm)]
    # Strings, for 'Stadnamn' interface: utm positions as a vector of strings like "3,233"
    vsutm = map(utm -> "$(utm[1]),$(utm[2])", vutm)
    sonline = get_config_value("Behaviour when data is missing", "Allow online geographical name collection", String)
    vfoundnames = point_names(vsutm; online = sonline == "false" ? false : true)
    vname = [n == "" ? sutm : n   for (n, sutm) in zip(vfoundnames, vsutm)]
    # Read 'surface bending' σ from csv files. This is for manual filtering or other adjustment of pararameters.
    vσ = read_σ_from_csvs(smb, vutm)
    vsheet_ij, vcell_ij, vσ, vname
end

function read_σ_from_csvs(smb, vutm)
    # Read from all existing csv files into a dictionary utm => σ (position => stress)
    # The .csv file may contain many more entries than the length of vutm
    dic = Dict{CartesianIndex, Float32}()
    for sb in smb
        fofo = full_folder_path(sb)
        ffna_sum = joinpath(fofo, SUMMITS_FNAM)
        vutm_f = CartesianIndex.(read_indices_from_column(ffna_sum, 3))
        vσ_f = Float64.(readdlm(ffna_sum, '\t')[2:end, 5])
        pairs = vutm_f .=> vσ_f
        foreach(pairs) do p
            push!(dic, p)
        end
    end
    #
    Float32.(map(vutm) do utm
        get(dic, CartesianIndex(utm), NaN32)
    end)
end