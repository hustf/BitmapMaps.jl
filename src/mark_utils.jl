# NOTE:
# This is work in progress. Also related: t_prominence and mark_utils.jl. We are trying to identify and mark 
# peaks from their prominence ("primary factor"), and while to reusing some functions from 'contour.jl'.
# Marking may be delegated to the svg.

# This places simple markers
function mark_square!(bwbuf::Matrix{Gray{Bool}}, center_idx::Int, side::Int)
    rows, cols = size(bwbuf)
    center = linear_to_cartesian(center_idx, cols)
    half_side = div(side, 2)

    # Define the square marking function
    function mark!(window)
        window .= Gray{Bool}(true)
    end

    # Create an offset array centered around the target pixel
    offset_bwbuf = OffsetArray(bwbuf, -half_side:half_side, -half_side:half_side)
    
    # Apply the marking function to the window around the center
    mapwindow!(mark!, offset_bwbuf, (half_side, half_side), center)
end

function mark_equilateral_triangle!(bwbuf::Matrix{Gray{Bool}}, base_center_idx::Int, side::Int)
    rows, cols = size(bwbuf)
    base_center = linear_to_cartesian(base_center_idx, cols)
    height = round(Int, side * sqrt(3) / 2)
    half_side = div(side, 2)

    # Define the triangle marking function
    function mark!(window, center)
        for i in 0:height
            row = center[1] - i
            if row < 1 || row > size(window, 1)
                continue
            end
            start_col = center[2] - div(i * side, height)
            end_col = center[2] + div(i * side, height)
            for col in start_col:end_col
                if col >= 1 && col <= size(window, 2)
                    window[row, col] = Gray{Bool}(true)
                end
            end
        end
    end

    # Create an offset array centered around the target pixel
    offset_bwbuf = OffsetArray(bwbuf, -height:0, -half_side:half_side)
    
    # Apply the marking function to the window around the base center
    mapwindow!((window) -> mark!(window, base_center), offset_bwbuf, (height, half_side), base_center)
end

# LinearIndices(A::AbstractArray)
# Return a LinearIndices array with the same shape and axes as A, holding the linear index of each entry in A. Indexing this array with cartesian indices allows mapping them to linear indices.
function linear_to_cartesian(index, cols)
    r, c = divrem(index - 1, cols)
    CartesianIndex(r + 1, c + 1)
end