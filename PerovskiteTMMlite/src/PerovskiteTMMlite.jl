module PerovskiteTMMlite

# =============================================================================
# PerovskiteTMMlite.jl
# -----------------------------------------------------------------------------
# A small, dependency-light transfer-matrix-method (TMM) optics package for
# multilayer thin-film solar cells (perovskite / perovskite-silicon tandems).
#
# It is a Julia translation of the coherent-TMM core of Steven Byrnes' Python
# `tmm` library (https://github.com/sbyrnes321/tmm, arXiv:1603.02720), plus a
# materials layer (n,k CSV handling) and a generation layer that turns the TMM
# field solution into a spatially- and spectrally-resolved carrier-generation
# rate G(x, λ).
#
# Public entry points (see each function's docstring for details):
#
#   Optics
#     coh_tmm(pol, n_list, d_list, th_0, λ)   -> Dict with R, T, absorption, fields
#     absorp_in_each_layer(coh_data)          -> per-layer absorptance vector
#     position_resolved(layer, z, coh_data)   -> Poynting / absorbed density / E
#
#   Materials
#     load_nk_data(path)                      -> NKData   (parse a wl,n,k CSV, nm)
#     interpolate_nk(nk, λ_grid)              -> Vector{ComplexF64}
#     load_and_interp(path, λ_grid)           -> Vector{ComplexF64}
#
#   Stacks (data/nk/*.csv must be present)
#     build_stack1()      air/MgF2/ITO/SnO2/C60/CsFAPbIBr/ITO_inter/Si   (paper optic)
#     build_stack2()      ... + spiro-TTB HTL, Si-terminated
#     build_stack1_Ag()   single-junction, Ag back reflector
#
#   Generation
#     compute_G_matrix(stack, λ_grid, x_grid_nm)   -> (G, A_per_layer, R, T)
#     tmm_generation_slice(λ, stack, grid_m; ...)  -> G on device nodes (SI)
#     beer_lambert_generation(...)                 -> Beer–Lambert reference
#
# Include order matters (Julia resolves method definitions per file at include
# time): tmm_core -> materials -> generation -> stacks. Do not reorder.
# =============================================================================

using LinearAlgebra        # matrix multiply used by coh_tmm (make it explicit,
                           #   so `using PerovskiteTMMlite` alone is self-sufficient)

include("tmm_core.jl")     # Fresnel + coh_tmm + list_snell + R_from_r + EPSILON
include("materials.jl")    # NKData, load_nk_data, interpolate_nk, load_and_interp
include("generation.jl")   # AbsorpAnalyticFn, absorp_in_each_layer, Stack,
                           #   compute_G_matrix, tmm_generation_slice, ...
include("stacks.jl")       # build_stack1 / build_stack2 / build_stack1_Ag

# ----------------------------- public API ------------------------------------
# Optics
export coh_tmm, absorp_in_each_layer, position_resolved
export interface_r, interface_t          # Fresnel primitives (handy for tests)

# Materials
export NKData, load_nk_data, interpolate_nk, load_and_interp

# Stacks
export Stack, build_stack1, build_stack2, build_stack1_Ag
export layer_nk_at, layer_starts_nm

# Generation
export compute_G_matrix, tmm_generation_slice, beer_lambert_generation
export AbsorpAnalyticFn, fill_in!, evaluate

end # module PerovskiteTMMlite
