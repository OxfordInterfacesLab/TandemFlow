module CTUtils

## Helper functions for simulations
export generation_from_scaps
export save_cell_profile, save_iv
export tanh_diffusion_profile
export get_cell_characteristics

## Functions for comparing SCAPS calculations and ChargeTransport calculations
# TODO: Figure out where these helpers are supposed to go
export IV
export CellProfile
export parse_scaps
export ct_to_profile
export set_IV!
export compare_densities, compare_bands, compare_iv

include("utils.jl")

end