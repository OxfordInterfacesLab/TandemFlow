#=

Script for simulation of a tandem solar cell

TOP CELL: PSC_C60_PVK_PTAA
BOTTOM CELL: SI_TOPCON
ARCHITECTURE: 4 TERMINAL

=#

include("../cells/Si_TOPCON.jl")
include("../cells/PSC_C60_PVK_PTAA.jl")
using .Si_TOPCON
using .PSC_2