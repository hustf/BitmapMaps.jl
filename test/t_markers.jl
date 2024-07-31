# Markers is a lightweight utility in BitmapMaps used for drawing on
# existing images.
# There's probably a library somewhere for this, too. Otherwise, 
# this might be separated out as it's own package.
# Note, the interface would naturally be better using symbols (:cross, etc...)
using Test
using BitmapMaps
using BitmapMaps: Gray, mark_at!
w, h = 100, 100
bw = zeros(Gray{Bool}, h, w)
x1, y1 = 5, 5
x2, y2 = 20, 50
# Pixel indexes - name hints at the y, then x sequence (row, column)
pis = [CartesianIndex((y1, x1)), CartesianIndex((y2, x2))]
p_nw = CartesianIndex((1,1))
########################
# Default: Hollow square
########################

mark_at!(bw, pis)
@test sum(bw) == 16 # 2 * (2 * 3 + 2)
# 
fill!(bw, false)
mark_at!(bw, pis, side = 9)
@test sum(bw) == 64 # 2 * (2 * 9 + 2 * (9 - 2)) 
#
fill!(bw, false)
mark_at!(bw, CartesianIndex(h ÷ 2, w ÷ 2), side = w - 1)
@test sum(bw) == 392 # 1 * (2 * 99 + 2 * (99 - 2)) 
# Over the edge
fill!(bw, false)
mark_at!(bw, p_nw)
@test sum(bw) == 3 


###############
# Hollow circle
###############

f_is_filled = BitmapMaps.func_is_on_circle
fill!(bw, false)
mark_at!(bw, pis; f_is_filled)
@test sum(bw) == 18
# 
fill!(bw, false)
mark_at!(bw, pis, side = 9; f_is_filled)
@test sum(bw) == 64 
#
fill!(bw, false)
mark_at!(bw, CartesianIndex(h ÷ 2, w ÷ 2), side = w - 1; f_is_filled)
@test sum(bw) == 312 
# Over the edge
fill!(bw, false)
mark_at!(bw, p_nw; f_is_filled)
@test sum(bw) == 4

#############################
# Hollow equilateral triangle
#############################

f_is_filled = BitmapMaps.func_is_on_triangle
fill!(bw, false)
mark_at!(bw, pis; f_is_filled)
@test sum(bw) == 6
# 
fill!(bw, false)
mark_at!(bw, pis, side = 9; f_is_filled)
@test sum(bw) == 46
#
fill!(bw, false)
mark_at!(bw, CartesianIndex(h ÷ 2, w ÷ 2), side = w - 1; f_is_filled)
@test sum(bw) == 313
# Over the edge
fill!(bw, false)
mark_at!(bw, p_nw; f_is_filled)
@test sum(bw) == 1

#################
# Horizontal line
#################

f_is_filled = BitmapMaps.func_is_on_hline
fill!(bw, false)
mark_at!(bw, pis; f_is_filled)
@test sum(bw) == 6
# 
fill!(bw, false)
mark_at!(bw, pis, side = 9; f_is_filled)
@test sum(bw) == 18
#
fill!(bw, false)
mark_at!(bw, CartesianIndex(h ÷ 2, w ÷ 2), side = w - 1; f_is_filled)
display_if_vscode(bw)
@test sum(bw) == 99

#################
# Vertical line
#################

f_is_filled = BitmapMaps.func_is_on_vline
fill!(bw, false)
mark_at!(bw, pis; f_is_filled)
@test sum(bw) == 6
# 
fill!(bw, false)
mark_at!(bw, pis, side = 9; f_is_filled)
@test sum(bw) == 18
#
fill!(bw, false)
mark_at!(bw, CartesianIndex(h ÷ 2, w ÷ 2), side = w - 1; f_is_filled)
display_if_vscode(bw)
@test sum(bw) == 99

#######
# Cross
#######

f_is_filled = BitmapMaps.func_is_on_cross
fill!(bw, false)
mark_at!(bw, pis; f_is_filled)
@test sum(bw) == 10 #  2 * (2 * 3 - 1)
# 
fill!(bw, false)
mark_at!(bw, pis, side = 9; f_is_filled)
@test sum(bw) == 34 # 2 * (2 * side - 1)
#
fill!(bw, false)
mark_at!(bw, CartesianIndex(h ÷ 2, w ÷ 2), side = w - 1; f_is_filled)
display_if_vscode(bw)
@test sum(bw) == 197 #  (2 * 99 - 1)

#########
# X-Cross
#########

f_is_filled = BitmapMaps.func_is_on_xcross
fill!(bw, false)
mark_at!(bw, pis; f_is_filled)
@test sum(bw) == 10
# 
fill!(bw, false)
mark_at!(bw, pis, side = 9; f_is_filled)
@test sum(bw) == 34
#
fill!(bw, false)
mark_at!(bw, CartesianIndex(h ÷ 2, w ÷ 2), side = w - 1; f_is_filled)
display_if_vscode(bw)
@test sum(bw) == 197

######################
# Equilateral triangle
######################

f_is_filled = BitmapMaps.func_is_in_triangle
fill!(bw, false)
mark_at!(bw, pis; f_is_filled)
@test sum(bw) == 8
# 
fill!(bw, false)
mark_at!(bw, pis, side = 9; f_is_filled)
@test sum(bw) == 50
#
fill!(bw, false)
mark_at!(bw, CartesianIndex(h ÷ 2, w ÷ 2), side = w - 1; f_is_filled)
display_if_vscode(bw)
@test sum(bw) == 3110

######
# Disc
######

f_is_filled = BitmapMaps.func_is_in_circle
fill!(bw, false)
mark_at!(bw, pis; f_is_filled)
@test sum(bw) == 18
# 
fill!(bw, false)
mark_at!(bw, pis, side = 9; f_is_filled)
@test sum(bw) == 138
#
fill!(bw, false)
mark_at!(bw, CartesianIndex(h ÷ 2, w ÷ 2), side = w - 1; f_is_filled)
display_if_vscode(bw)
@test sum(bw) == 7705

########
# Square
########

f_is_filled = BitmapMaps.func_is_in_square
fill!(bw, false)
mark_at!(bw, pis; f_is_filled)
@test sum(bw) == 18
# 
fill!(bw, false)
mark_at!(bw, pis, side = 9; f_is_filled)
@test sum(bw) == 162
#
fill!(bw, false)
mark_at!(bw, CartesianIndex(h ÷ 2, w ÷ 2), side = w - 1; f_is_filled)
display_if_vscode(bw)
@test sum(bw) == 99^2

#################
# String argument
#################
fill!(bw, false)
mark_at!(bw, pis, 3, "in_square")
@test sum(bw) == 18