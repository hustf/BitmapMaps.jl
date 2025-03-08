# Used by steps in the pipeline.
#
# Utilties concerning ImageSegmentation.SegmentedImage.
# Functions relating segmentation of an image with data,
# primarily dictionaries of those data.

function segment_dictionary(M, segments::SegmentedImage, operator)
    @assert axes(labels_map(segments)) == axes(M)
    dic = Dict{Int64, eltype(M)}()
    R = CartesianIndices(axes(labels_map(segments)))
    for ci in R
        label = labels_map(segments)[ci]
        previous_value = get(dic, label, nothing)
        if isnothing(previous_value)
            value = M[ci]
        else
            value = operator(previous_value, M[ci])
        end
        push!(dic, label => value)
    end
    dic
end

function cumulative_cartesian_indices_dictionary(segments)
    R = CartesianIndices(axes(labels_map(segments)))
    segment_dictionary(R, segments, +)
end

function mean_index_dictionary(segments)
    R = CartesianIndices(axes(labels_map(segments)))
    # All dicts indexed by segment no.
    diccum = cumulative_cartesian_indices_dictionary(segments)
    dicmean = Dict{Int64, CartesianIndex}()
    for label in segment_labels(segments)
        cum = diccum[label]
        count = segment_pixel_count(segments, label)
        mean = CartesianIndex{2}((cum.I .÷ count)) # We stick to integer indices here.
        # Checking for error
        if mean ∉ R
            @show mean cum count R
            throw("Find out why")
        end
        push!(dicmean, label => mean)
    end
    dicmean
end

function bbox_internal_indices_dictionary(segments)
    R = CartesianIndices(axes(labels_map(segments)))
    # Top left, bottom right corner indices
    tld = segment_dictionary(R, segments, min)
    brd = segment_dictionary(R, segments, max)
    dic = Dict{Int64, CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}}()
    for label in segment_labels(segments)
        tl = tld[label]
        br = brd[label]
        bx = CartesianIndices((tl[1]:br[1], tl[2]:br[2]))
        push!(dic, label => bx)
    end
    dic
end

bbox_diagonal_length_dictionary(segments, bbdic) = Dict([label=> hypot((bb[end] - bb[1]).I...) for (label, bb) in bbdic])
