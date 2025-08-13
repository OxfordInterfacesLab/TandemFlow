module CTUtils

## Helper functions for simulations
export generation_from_scaps, scaps_to_df_generation
export save_iv

## Functions for comparing SCAPS calculations and ChargeTransport calculations
export CellInfo
export scaps_to_df, scaps_to_df_iv
export ct_to_df
export compare_densities, compare_bands, compare_iv
export find_info_scaps, find_info_ct

include("scaps_utils.jl")
include("sim_utils.jl")

end