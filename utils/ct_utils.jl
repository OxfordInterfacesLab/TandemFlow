module CTUtils


## Functions for comparing SCAPS calculations and ChargeTransport calculations
export CellInfo
export scaps_to_df
export ct_to_df
export compare_densities
export compare_bands
export find_info_scaps, find_info_ct

include("scaps_utils.jl")

end